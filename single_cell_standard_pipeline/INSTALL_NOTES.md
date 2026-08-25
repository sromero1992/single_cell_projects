# INSTALL NOTES — the packages that actually break

`00_rlibs_installation.R` handles everything on CRAN, Bioconductor, and the well-behaved GitHub packages. This file covers the ones that reliably cause trouble, in rough order of how much pain they cause.

**Order matters.** Do sections 1–4 in R. Section 5 (Python) is only needed for `10_trajectory_cellrank.R` — nothing else in the pipeline requires Python.

---

## Quick reference

| Package | Needed by | Install route | Difficulty |
|---|---|---|---|
| DoubletFinder | Script 01 (default) | GitHub | Easy |
| BPCells | Script 01 (optional) | GitHub, needs compiler | Medium |
| SCEVAN | Script 01 (optional) | GitHub | Medium |
| **CytoTRACE 2** | **Script 09** | **GitHub, subdir** | **Medium** |
| **CytoTRACE v1** | **Script 09 (optional)** | **GitHub, stale deps** | **Hard** |
| CellChat | Script 08 | GitHub | Medium |
| monocle3 | optional | GitHub, heavy system deps | Hard |
| **scanpy + CellRank** | **Script 10** | **conda/pip** | **Medium** |

---

## 0. Before anything else

```r
install.packages(c("devtools", "remotes", "BiocManager"))
options(timeout = 1200)   # GitHub builds routinely exceed the 60s default
```

The `timeout` line matters more than it looks — the single most common "install failed" report is a silent download timeout on a large tarball. `00_rlibs_installation.R` already sets it.

On a Linux/HPC box, confirm you have a working toolchain before starting:

```bash
gcc --version && g++ --version && make --version
```

If any of those are missing, every source install below will fail with confusing errors about `.so` files.

### The `TamuScDSC` package (preprocessing engine for Scripts 01 / 01b / 02)

Separate from the packages below (which are external and prone to breaking),
the pipeline's own preprocessing code is shipped as a local package, `TamuScDSC/`.
`00_rlibs_installation.R` installs it automatically as its **last step** (section
6), after its dependencies. You only need to do this by hand if the auto-install
couldn't locate the folder:

```r
devtools::install("TamuScDSC")     # from the repo root (folder containing TamuScDSC/)
```

It builds in seconds and its heavy dependencies are only *Suggests* (already
covered by `00_rlibs_installation.R`). Full details + the three usage scenarios:
`TamuScDSC/README.md` and `TamuScDSC/ARCHITECTURE.md`, or the "Installing the TamuScDSC
package" section of `QUICK_START.md`.

---

## 1. CytoTRACE 2 — required by Script 09

The primary potency method. Pure R, **no Python required**.

```r
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(CytoTRACE2)
```

The `subdir` argument is mandatory — the R package lives in a subdirectory of the repo, and omitting it produces a confusing "does not appear to be a package" error.

Verify:

```r
packageVersion("CytoTRACE2")   # expect 1.1.0 or newer
```

For Script 09's Ensembl → symbol step (only triggered if your object carries
Ensembl rownames) also install the annotation packages — Script 00 does this
automatically:

```r
BiocManager::install(c("AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db"))
```

### Exact working call (validated) — how Script 09 invokes CytoTRACE 2

These arguments were validated on the scDVEP benchmark (Nestorowa 2016 HSC data,
1920 cells × 20475 genes) and are the ones Script 09 now uses. Each item marked
below caused failed runs when changed, so do not "simplify" them:

```r
## Input MUST be a dense data.frame (genes × cells) with gene-SYMBOL rownames.
## Sparse matrices cause SILENT failures; Ensembl rownames collapse every score.
if (!is.data.frame(ct2_input))
  ct2_input <- as.data.frame(as.matrix(ct2_input))

ct2_result <- CytoTRACE2::cytotrace2(
  ct2_input,
  species               = "mouse",   # or "human"
  is_seurat             = FALSE,      # pass the data.frame, NOT the Seurat object
  parallelize_models    = FALSE,      # avoids parallel-backend issues
  parallelize_smoothing = FALSE,      # REQUIRED: avoids socket-cluster hangs
  seed                  = 42
)
```

Pitfalls that repeatedly broke runs:

- **Do NOT pass `verbose = TRUE`** — it is not a valid argument and throws
  `unused argument (verbose = TRUE)`.
- **Keep `parallelize_smoothing = FALSE`.** `TRUE` reintroduces socket-cluster
  hangs (required on Windows, safe and recommended on Linux).
- **Dense data.frame only.** A sparse `dgCMatrix` produces no error but returns
  meaningless scores.
- **Gene symbols, not Ensembl IDs.** Script 09 maps Ensembl → symbol via
  `org.Mm.eg.db` / `org.Hs.eg.db` before calling `cytotrace2()`, then
  deduplicates symbols (keeping the highest mean-expression row).

Return value is a `data.frame` (rows = cells): `CytoTRACE2_Score`,
`CytoTRACE2_Potency`, `CytoTRACE2_Relative`, `preKNN_CytoTRACE2_Score`,
`preKNN_CytoTRACE2_Potency`. The continuous score is clipped to `[0,1]`.

A self-contained reference implementation ships as `run_cytotrace2.R` in the
project root (adds AUCell cell-cycle scoring and optional CellRank trajectories).

### Dependencies it pulls in

`data.table`, `doParallel`, `dplyr`, `ggplot2`, `HiClimR`, `magrittr`, `Matrix`, `plyr`, `Rfast`, `RSpectra`, `Seurat`, `SeuratObject`, `stringr`.

`Rfast` and `RSpectra` compile from source and need a C++ toolchain. `HiClimR` needs `ncdf4`, which needs the system NetCDF library:

```bash
# RHEL / Rocky / CentOS
sudo dnf install netcdf-devel gsl-devel
# Debian / Ubuntu
sudo apt install libnetcdf-dev libgsl-dev
```

### Known failure modes

| Symptom | Cause | Fix |
|---|---|---|
| `does not appear to be a package` | `subdir` omitted | Add `subdir = "cytotrace2_r"` |
| Seurat v4 + Matrix 1.6 conflict | Documented incompatibility | Upgrade Seurat, or downgrade Matrix to 1.5-4.1 |
| `Rfast` fails to compile | No C++17 toolchain | Install `gcc-c++`, retry |
| Killed / OOM during prediction | Too many parallel workers | Lower `CT2_NCORES` to 1–2 in Script 09 |
| `is.atomic(y) is not TRUE` | Input not a plain dense data.frame | Script 09 passes a dense data.frame with `is_seurat = FALSE`; check the input isn't sparse |
| Scores all identical / meaningless | Sparse input, or Ensembl rownames | Use a dense data.frame with gene-symbol rownames (Script 09 handles both) |
| `unused argument (verbose = TRUE)` | `verbose` is not a valid arg | Remove it |
| Hangs during smoothing | `parallelize_smoothing = TRUE` | Set it to `FALSE` |

### `'configure' exists but is not executable` (ncdf4 → HiClimR → CytoTRACE2)

CytoTRACE 2 depends on `HiClimR`, which depends on `ncdf4`, which builds from source. On HPC systems `ncdf4` commonly fails with:

```
ERROR: 'configure' exists but is not executable
   -- see the 'R Installation and Administration Manual'
```

The package is fine. The `configure` script lost its executable bit when the tarball was unpacked, which happens when `TMPDIR` sits on a filesystem mounted **`noexec`**, or on NFS that strips permission bits. A tell-tale sign in the same log:

```
Warning: invalid uid value replaced by that for user 'nobody'
```

That is NFS ownership squashing — strong evidence `TMPDIR` is on a network home directory.

> ⚠️ `00_rlibs_installation.R` sets `TMPDIR` to `$HOME/Rtemp`. If `$HOME` is NFS-mounted, that setting *causes* this failure. Override it for the install.

**First, find out where R is actually building.** This matters because two common assumptions are wrong:

```r
tempdir()                       # where the build REALLY happens
Sys.getenv("TMPDIR")
```

> ⚠️ **`Sys.setenv(TMPDIR = ...)` inside a running R session does NOT move the build directory.** R fixes `tempdir()` once, at startup. Changing `TMPDIR` afterwards has no effect on the current session — you will see the build still running under the original `/tmp/Rtmp<xxxx>` path. To change it, set `TMPDIR` in `~/.Renviron` (or export it in the shell) and **restart R**.

**Step 1 — is the build directory mounted `noexec`?** This is the usual cause, and hardened HPC systems very often mount `/tmp` that way:

```bash
findmnt -no TARGET,OPTIONS -T /tmp
findmnt -no TARGET,OPTIONS -T $HOME
findmnt -no TARGET,OPTIONS -T /var/tmp
```

Any of them containing `noexec` cannot run `configure`, no matter what the file permissions say.

**Step 2 — point R at a directory that permits execution.** Pick one from above that is *not* `noexec`, then set it before R starts. In `~/.Renviron`:

```
TMPDIR=/var/tmp
```

Restart R, confirm with `tempdir()`, then `install.packages("ncdf4")`.

**Step 3 — if every candidate is `noexec`, build somewhere you control:**

```bash
mkdir -p ~/rbuild && cd ~/rbuild
wget https://cran.r-project.org/src/contrib/ncdf4_1.24.tar.gz
tar xzf ncdf4_1.24.tar.gz
chmod +x ncdf4/configure
TMPDIR=~/rbuild R CMD INSTALL ncdf4
```

**Fix — unpack and restore the bit manually:**

```bash
cd /tmp
wget https://cran.r-project.org/src/contrib/ncdf4_1.24.tar.gz
tar xzf ncdf4_1.24.tar.gz
chmod +x ncdf4/configure
R CMD INSTALL ncdf4
```

> 💡 **`/tmp` mounted `noexec` breaks every source package with a `configure` script**, not just `ncdf4`. If `findmnt` shows `noexec` on `/tmp`, fix it once and permanently by putting an exec-capable directory in `~/.Renviron`:
> ```
> TMPDIR=/var/tmp
> ```
> Restart R and confirm with `tempdir()`. This prevents the whole class of failure.

### `nc-config not found` (the next error after the exec fix)

Once `configure` actually runs, it fails differently if the NetCDF library itself is absent:

```
Error, nc-config not found or not executable.
```

`nc-config` ships with the NetCDF development package. Check what you already have before installing anything:

```bash
which nc-config
rpm -q netcdf netcdf-devel     # RHEL / Rocky
module avail netcdf            # HPC module systems often provide it
```

**Option A — system package (needs root):**

```bash
sudo dnf install netcdf-devel        # RHEL / Rocky
sudo apt install libnetcdf-dev       # Debian / Ubuntu
```

**Option B — HPC module (no root needed):**

```bash
module load netcdf
which nc-config                      # confirm it is now on PATH
R CMD INSTALL ncdf4
```

**Option C — conda (no root needed).** Works, but read the caveat:

```bash
conda install -n scanpy_env_311 -c conda-forge libnetcdf
NC=$HOME/miniconda3/envs/scanpy_env_311/bin/nc-config
cd ~/rbuild
R CMD INSTALL --configure-args="--with-nc-config=$NC" ncdf4
```

> ⚠️ **Caveat:** linking an R package against conda's NetCDF means the compiled `ncdf4.so` needs conda's libraries at **run time**, reintroducing the same library-path fragility as the llvmlite/libstdc++ issue in section 5. If `library(ncdf4)` later fails with a missing `.so`, that is why. Prefer Options A or B when available.

**Option D — skip it.** This whole chain exists only for CytoTRACE 2. Set `RUN_CYTOTRACE2 <- FALSE` in Script 09; the entropy metrics have no dependencies and still produce a potency readout. This is a perfectly reasonable place to stop.

Then continue up the chain:

```r
install.packages("HiClimR")
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
```

**Fix 3 — conda-forge binaries, no compilation at all:**

```bash
conda install -c conda-forge r-ncdf4 r-hiclimr
```

Only useful if your R *is* the conda R. Check with `R.home()`; if it points inside `miniconda3`, this works.

**If none of it works:** set `RUN_CYTOTRACE2 <- FALSE` in Script 09. The entropy metrics have no dependencies and still run, so you keep a potency readout.

> **Note on the update prompt.** When `install_github` asks which packages to update, answering `3: None` (as is usually right) is fine — the `ncdf4` failure is unrelated to those updates.

### Conda alternative

If dependency resolution is fighting you, the authors ship an environment file:

```bash
git clone https://github.com/digitalcytometry/cytotrace2.git
cd cytotrace2
conda env create -f environment_R.yml
conda activate cytotrace2
R -e 'devtools::install_local("./cytotrace2_r")'
```

Takes 5–10 minutes normally, occasionally up to an hour. On M1/M2 Macs the conda solve fails; use:

```bash
conda create -n cytotrace2 && conda activate cytotrace2
conda config --env --set subdir osx-64
conda env update --file environment_R.yml
```

---

## 2. CytoTRACE v1 — optional, Script 09

The original method, kept as an independent cross-check. **This is the hardest install in the pipeline** — the package is from 2020 and its dependencies have moved on. Script 09 runs fine without it (`RUN_CYTOTRACE1 <- FALSE`), so do not lose a day to this one.

```r
devtools::install_github("gunsagargulati/CytoTRACE")
```

### Dependencies that break

```r
# sva is Bioconductor, not CRAN — install first or the build fails
BiocManager::install("sva")
install.packages(c("nnls", "egg", "ccaPP", "HiClimR", "plyr"))
```

`ccaPP` needs `RcppArmadillo`, which needs a BLAS/LAPACK dev package:

```bash
sudo dnf install openblas-devel lapack-devel     # RHEL family
sudo apt install libopenblas-dev liblapack-dev   # Debian family
```

### The Python part — only for `iCytoTRACE()`

The `iCytoTRACE()` function (multi-batch integration) needs Python with `scanoramaCT` and `numpy`. **Script 09 does not call `iCytoTRACE()`** — it runs the plain `CytoTRACE()` per sample instead, which avoids this entirely. Skip this unless you specifically want integrated CytoTRACE:

```bash
pip install numpy scanoramaCT
```

### Known failure modes

| Symptom | Fix |
|---|---|
| `there is no package called 'sva'` | `BiocManager::install("sva")` — it's Bioconductor, not CRAN |
| Fails on modern Matrix | Known; the package predates Matrix 1.5+. Use `RUN_CYTOTRACE1 <- FALSE` |
| Enormous memory use | `CytoTRACE()` densifies the matrix. Script 09 already runs it per sample; reduce `MAX_CELLS` if needed |

If it will not install after 30 minutes of effort, set `RUN_CYTOTRACE1 <- FALSE`. CytoTRACE 2 plus the built-in entropy metrics already give you two independent estimates.

---

## 2b. CCAT / SCENT — optional potency benchmark, Script 09

CCAT (Teschendorff signalling-entropy connectome correlation) is an optional
potency score in Script 09 (`RUN_CCAT`). SCENT's slower signalling-entropy rate
is behind `RUN_SCENT` (off by default). Both are **fully guarded** — if any part
of the install or the gene-ID chain fails, Script 09 prints `[SKIP]` and finishes
without them.

**Install (Script 00 does all of this):**

```r
remotes::install_version("smoother", version = "1.1")  # archived on CRAN — pin it
devtools::install_github("aet21/SCENT")                # provides net17Jan16.m, CompCCAT, DoIntegPPI
install.packages(c("homologene", "mclust", "pROC"))    # ortholog map + benchmark helpers
BiocManager::install(c("org.Hs.eg.db", "biomaRt", "scuttle"))
```

`smoother` must go in **before** SCENT (SCENT depends on it), which is why
Script 00 installs the pinned version first.

**The gene-ID chain (this is what breaks if skipped).** CCAT correlates each
cell's expression with PPI hub degree, and SCENT's PPI (`net17Jan16.m`) is indexed
by **human Entrez**. So mouse data must be mapped:

```
mouse symbol  --homologene(inTax=10090, outTax=9606)-->  human symbol
              --org.Hs.eg.db-->                          human Entrez
              --match into net17Jan16.m-->               PPI hub degree
```

Script 09's STEP 3c does exactly this (`CCAT_SPECIES = "mouse"`), then
`DoIntegPPI` → `CompCCAT`. Set `CCAT_SPECIES = "human"` to skip the ortholog hop.

**Known failure modes:**

| Symptom | Cause | Fix |
|---|---|---|
| `smoother` won't install | archived on CRAN | `remotes::install_version("smoother","1.1")` |
| CCAT drops all genes | mouse symbols vs human PPI | ensure `homologene` installed; `CCAT_SPECIES="mouse"` |
| `[SKIP] Need SCENT + org.Hs.eg.db` | a dependency is missing | install the four blocks above |
| CCAT inverts on quiescent cells | it's a breadth method (HSC paradox) | expected biology — the concordance step flags it |
| SCENT (`CompSR`) very slow | signalling-entropy is O(min) on 1000+ cells | leave `RUN_SCENT = FALSE` unless you need it |

CCAT/SCENT are treated as potency methods: their `CCAT_score` / `SCENT_score`
join the Spearman concordance, the by-condition plots, and the per-contrast
comparisons in Script 09.

---

## 3. DoubletFinder, BPCells, SCEVAN — Script 01

### DoubletFinder (default doublet caller)

```r
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

Usually painless. Note the function names changed in the 2023 rewrite: `paramSweep_v3` → `paramSweep`, `doubletFinder_v3` → `doubletFinder`. **The unified pipeline uses the NEW names.** If you get `could not find function "paramSweep"`, you have an old version:

```r
remove.packages("DoubletFinder")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)
```

### BPCells (only if `USE_BPCELLS <- TRUE`)

```r
remotes::install_github("bnprks/BPCells/r")
```

Needs HDF5 development headers:

```bash
sudo dnf install hdf5-devel     # RHEL family
sudo apt install libhdf5-dev    # Debian family
```

If it will not build, set `USE_BPCELLS <- FALSE` in Script 01 — it is a memory optimisation, not a requirement.

### SCEVAN (only if `RUN_SCEVAN <- TRUE`)

```r
remotes::install_github("AntonioDeFalco/SCEVAN")
```

Pulls `yaGST` from GitHub. If `ggtree` fails, install it from Bioconductor instead of GitHub:

```r
BiocManager::install("ggtree")
```

Set `RUN_SCEVAN <- FALSE` if your study has no tumour samples — it is skippable.

---

## 4. CellChat and monocle3 — Scripts 08 and optional

### CellChat

```r
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors", "Biobase"))
remotes::install_github("jinworks/CellChat")
```

`ComplexHeatmap` must come from Bioconductor first; installing CellChat before it produces a misleading error about `circlize`.

### monocle3 (not used by any script; optional)

The heaviest install here. System libraries first:

```bash
# RHEL family
sudo dnf install gdal-devel proj-devel geos-devel udunits2-devel sqlite-devel
# Debian family
sudo apt install libgdal-dev libproj-dev libgeos-dev libudunits2-dev libsqlite3-dev
```

Then:

```r
BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "HDF5Array",
                       "terra", "ggrastr"))
remotes::install_github("cole-trapnell-lab/leidenbase")
remotes::install_github("cole-trapnell-lab/monocle3")
```

`sf` is the usual failure point and needs the GDAL config flags that `00_rlibs_installation.R` already passes.

---

## 5. Python environment for CellRank — Script 10 only

**Only `10_trajectory_cellrank.R` needs this.** Scripts 00–09 are pure R. If you are not doing trajectory inference, stop here.

CellRank and scanpy are Python-only. Script 10 drives them from R through `reticulate`.

### 5.0 CPU only — no GPU, no PyTorch

Nothing in this pipeline uses a GPU, and **PyTorch is not a dependency**:

| Component | Backend | GPU? |
|---|---|---|
| CellRank `CytoTRACEKernel` (Script 10) | pyGPCCA / scipy | No |
| scanpy / scvelo preprocessing | numpy / scipy | No |
| CytoTRACE 2 (Script 09) | **R package**, `.rds` model matrices | No |

Script 09's CytoTRACE 2 is the *R* package and needs no Python at all. Only Script 10 does.

### 5.1 Automated setup (recommended)

Use `setup_python_env.sh`, which installs into an existing conda env, snapshots it first for rollback, verifies every import, and prints the reticulate lines for R:

```bash
chmod +x setup_python_env.sh
./setup_python_env.sh --check      # verify only, install nothing
./setup_python_env.sh              # install into scanpy_env_311
./setup_python_env.sh --optional   # also install cellmapper
./setup_python_env.sh --env my_env # target a different environment
```

Edit `ENV_NAME` at the top of the script to change the default target. It is safe to re-run — pip skips anything already satisfied.

### 5.2 Manual equivalent

```bash
conda activate scanpy_env_311
pip install "cellrank>=2.0" scvelo
python -c "import cellrank, scanpy, scvelo; print(cellrank.__version__, scanpy.__version__)"
```

Optional, only if you set `CT_IMPUTE_METHOD <- "cellmapper"`:

```bash
pip install cellmapper
```

Verify in the shell before touching R — debugging import failures through the reticulate layer is considerably worse.

> **Adding to a working env carries some risk.** pip can pull a numpy/scipy version that conflicts with what scanpy already has. The setup script writes a `pip freeze` snapshot first so you can restore. To avoid the risk entirely, clone first:
> ```bash
> conda create -n scanpy_cellrank --clone scanpy_env_311
> ```

### 5.3 Building a fresh environment instead

```bash
conda create -n cellrank_env python=3.11 -y
conda activate cellrank_env
pip install "cellrank>=2.0" scanpy scvelo anndata igraph leidenalg
```

### 5.2 Install reticulate in R

```r
install.packages("reticulate")
```

### 5.3 Point Script 10 at the environment

In `10_trajectory_cellrank.R`, Section 1.1:

```r
PYTHON_ENV_TYPE <- "conda"
PYTHON_ENV      <- "cellrank_env"
```

Script 10 binds the interpreter and verifies every required module **before** doing any work, so a misconfiguration fails in seconds with a clear message rather than after loading a large object.

### 5.4 Sanity check from R

```r
library(reticulate)
use_condaenv("cellrank_env", required = TRUE)
py_config()                          # confirm the path is the env you built
py_module_available("cellrank")      # must be TRUE
```

### 5.4b `py_config()` looks right but `py_module_available()` is FALSE

This is the single most common conda + reticulate failure on RHEL/Rocky, and it looks nothing like a missing package:

```r
py_config()                        # correct interpreter, correct numpy
py_module_available("cellrank")    # FALSE
```

**`py_module_available()` returns FALSE when the import RAISES**, not only when the module is absent. It swallows the exception. The module is installed; importing it inside R is what fails.

**Step 1 — see the actual error. Do not skip this.**

```r
reticulate::py_run_string("import cellrank")
```

**Step 2 — read the message.**

| Error contains | Cause |
|---|---|
| `GLIBCXX_3.4.xx not found` / `libstdc++.so.6` | R and conda disagree on libstdc++ (most common) |
| `numpy.core.multiarray failed to import` | numpy ABI mismatch |
| `undefined symbol` | compiled-extension mismatch, same family |

**The confirmed real-world error looks like this:**

```
OSError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found
  (required by .../site-packages/llvmlite/binding/../../../../libLLVM-14.so)
```

Chain: `cellrank → scanpy → numba → llvmlite → libLLVM-14.so → libstdc++`. The system `/lib64/libstdc++.so.6` is older than what conda's `libLLVM` needs. Nothing is missing; the wrong copy of one library is being loaded.

**Step 3a — the surgical fix (preferred): `LD_PRELOAD` just libstdc++.**

Replacing only that one library is safer than reordering the whole path, because R itself was compiled against the system libraries:

```bash
export LD_PRELOAD=/home/ssromerogon/miniconda3/envs/scanpy_env_311/lib/libstdc++.so.6
R
```

Persist it in `~/.Renviron`:

```
LD_PRELOAD=/home/ssromerogon/miniconda3/envs/scanpy_env_311/lib/libstdc++.so.6
RETICULATE_PYTHON=/home/ssromerogon/miniconda3/envs/scanpy_env_311/bin/python
```

First confirm conda's copy actually has the required version:

```bash
strings /home/ssromerogon/miniconda3/envs/scanpy_env_311/lib/libstdc++.so.6 \
  | grep GLIBCXX_3.4.30
```

Empty output means conda's is also too old — install a newer one:

```bash
conda install -n scanpy_env_311 -c conda-forge libstdcxx-ng
```

**Step 3b — the broader fix: `LD_LIBRARY_PATH`.**

System R on RHEL loads the OS `libstdc++`, which is older than the one conda's compiled extensions (numba, llvmlite, scipy) were built against. The conda `lib` directory must come **first** on the library path.

> ⚠️ **Check this first if you have run `00_rlibs_installation.R` in the same session.** That script's "FORCE protocol" sets
> ```r
> Sys.setenv(LD_LIBRARY_PATH = "/usr/lib64")
> ```
> which pins the linker to system libraries — exactly the condition that breaks conda extensions. Verify with `Sys.getenv("LD_LIBRARY_PATH")`. If it reads `/usr/lib64`, that is very likely your problem. Use a **fresh R session** that has not sourced Script 00.

`LD_LIBRARY_PATH` is read by the dynamic linker when the **process starts**, so `Sys.setenv()` inside a running R session is unreliable. Set it before R launches. Put this in `~/.Renviron`:

```
LD_LIBRARY_PATH=/home/ssromerogon/miniconda3/envs/scanpy_env_311/lib:/usr/lib64
RETICULATE_PYTHON=/home/ssromerogon/miniconda3/envs/scanpy_env_311/bin/python
```

Then **fully restart R** (not just `use_condaenv` again — reticulate binds once per session).

Equivalent, from a shell:

```bash
export LD_LIBRARY_PATH=/home/ssromerogon/miniconda3/envs/scanpy_env_311/lib:$LD_LIBRARY_PATH
export RETICULATE_PYTHON=/home/ssromerogon/miniconda3/envs/scanpy_env_311/bin/python
R
```

**Step 4 — if that is not enough**, give conda a newer libstdc++:

```bash
conda install -n scanpy_env_311 -c conda-forge libstdcxx-ng
```

**Step 5 — last resort**, use a virtualenv instead of conda. Virtualenvs use the system Python and therefore the system libstdc++, so the mismatch cannot arise:

```bash
python3 -m venv ~/cellrank_venv
source ~/cellrank_venv/bin/activate
pip install "cellrank>=2.0" scanpy scvelo anndata igraph leidenalg
```

```r
PYTHON_ENV_TYPE <- "virtualenv"
PYTHON_ENV      <- "~/cellrank_venv"
```

**Confirming the fix:**

```r
library(reticulate)
use_condaenv("scanpy_env_311", required = TRUE)
cr <- import("cellrank")
cr$`__version__`                      # "2.0.7"
py_module_available("cellrank")       # TRUE
```

Note `import()` is the better test than `py_module_available()` — it shows you the error instead of hiding it behind FALSE.

### 5.5 Known failure modes

| Symptom | Cause | Fix |
|---|---|---|
| `Unable to locate conda environment` | reticulate can't find conda | Set `options(reticulate.conda_binary = "/path/to/conda")`, or use `PYTHON_ENV_TYPE <- "python"` with an explicit `PYTHON_BIN` |
| Python bound to the wrong interpreter | Something imported Python before `use_condaenv()` | **Restart R.** reticulate binds once per session and cannot be switched afterwards |
| `ImportError: numpy.core.multiarray failed` | NumPy 1.x / 2.x ABI mismatch | Rebuild the env cleanly; do not mix conda and pip for numpy |
| Segfault on import | conda R and conda Python fighting over libstdc++ | Use system R with a conda *Python*, not conda R |
| `CytoTRACEKernel` missing | CellRank v1 installed | `pip install --upgrade "cellrank>=2.0"` — the kernel API changed in v2 |
| `compute_cytotrace()` complains about the layer | No imputed layer present | Script 10 computes it via `CT_IMPUTE_METHOD`; keep it at `"moments"` |
| Cannot write `.h5ad` | Enum columns in `.obs` (known CellRank issue) | Non-fatal — the `.rds` still saves. Set `SAVE_H5AD <- FALSE` |

### 5.6 HPC note

On a cluster, load the module stack before creating the environment, and pin `reticulate` to the environment inside your job script:

```bash
module load gcc/11 hdf5
export RETICULATE_PYTHON=/home/$USER/miniconda3/envs/cellrank_env/bin/python
Rscript 10_trajectory_cellrank.R
```

`RETICULATE_PYTHON` overrides everything else and is the most reliable way to get a deterministic binding in a batch job.

---

## 6. Verifying the whole environment

Run this after finishing the sections you need. It mirrors the gate in `00_rlibs_installation.R` but also covers the Script 09/10 additions:

```r
check <- function(pkgs) {
  st <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  data.frame(Package = pkgs, Installed = st, row.names = NULL)
}

cat("--- Core (Scripts 01-08) ---\n")
print(check(c("Seurat", "SeuratWrappers", "harmony", "celda",
              "DoubletFinder", "scDblFinder", "SCEVAN",
              "MAST", "CellChat", "enrichR", "writexl", "openxlsx")))

cat("\n--- Script 09 (cell scores) ---\n")
print(check(c("CytoTRACE2", "CytoTRACE")))

cat("\n--- Script 10 (trajectory) ---\n")
print(check("reticulate"))
if (requireNamespace("reticulate", quietly = TRUE)) {
  try({
    reticulate::use_condaenv("cellrank_env", required = TRUE)
    for (m in c("scanpy", "cellrank", "scvelo", "anndata")) {
      cat(sprintf("  %-10s %s\n", m,
          if (reticulate::py_module_available(m)) "OK" else "MISSING"))
    }
  }, silent = TRUE)
}
```

### What is genuinely required vs. skippable

| Must have | Can skip |
|---|---|
| Seurat, dplyr, ggplot2, writexl | BPCells (`USE_BPCELLS <- FALSE`) |
| DoubletFinder **or** scDblFinder | SCEVAN (`RUN_SCEVAN <- FALSE`) |
| celda (DecontX) | CytoTRACE v1 (`RUN_CYTOTRACE1 <- FALSE`) |
| CytoTRACE2 — *for Script 09* | monocle3 (unused) |
| reticulate + Python stack — *for Script 10* | Script 10 entirely, if not doing trajectories |

The entropy metrics in Script 09 have **no external dependencies** and always run, so you get a potency readout even if every optional install above fails.

---

## Sources

- [CytoTRACE 2 — digitalcytometry/cytotrace2](https://github.com/digitalcytometry/cytotrace2)
- [CytoTRACE 2 R README](https://github.com/digitalcytometry/cytotrace2/blob/main/README.md)
- [Kang et al., *Nature Methods* 2025 — CytoTRACE 2](https://www.nature.com/articles/s41592-025-02857-2)
- [CytoTRACE v1 — gunsagargulati/CytoTRACE](https://github.com/gunsagargulati/CytoTRACE)
- [CellRank — CytoTRACEKernel API](https://cellrank.readthedocs.io/en/latest/api/_autosummary/kernels/cellrank.kernels.CytoTRACEKernel.html)
- [CellRank meets CytoTRACE tutorial](https://cellrank.readthedocs.io/en/stable/notebooks/tutorials/kernels/400_cytotrace.html)
