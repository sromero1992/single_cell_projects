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
| `is.atomic(y) is not TRUE` | Input not a plain matrix/Seurat | Script 09 passes a Seurat object with `is_seurat = TRUE`; check your object isn't wrapped |

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

### 5.1 Create the environment

Do this **once**, from a terminal — not from inside R, and not from inside an analysis script:

```bash
# mamba is much faster than conda for this solve; either works
conda create -n cellrank_env python=3.11 -y
conda activate cellrank_env

pip install "cellrank>=2.0" scanpy scvelo anndata igraph leidenalg
```

Optional, for the newer imputation route (`CT_IMPUTE_METHOD <- "cellmapper"`):

```bash
pip install cellmapper
```

Verify before touching R:

```bash
python -c "import cellrank, scanpy, scvelo; print(cellrank.__version__, scanpy.__version__)"
```

If that line errors, fix it here. Debugging Python import failures through the reticulate layer is considerably worse.

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
