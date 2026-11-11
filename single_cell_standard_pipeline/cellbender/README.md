# CellBender ambient-RNA removal (per sample)

Optional preprocessing step that removes ambient RNA / empty-droplet background
from each sample's **raw** 10x matrix with
[CellBender](https://github.com/broadinstitute/CellBender) `remove-background`,
**before** the R pipeline (Script 01). It is an alternative to the built-in
**DecontX** step — use one or the other, not both.

## How it fits the pipeline (short answer to the usual questions)

- **It runs externally, in a Linux terminal**, in its own conda env (GPU strongly
  recommended). CellBender is a PyTorch program; it is **not** part of the R/
  TamuScDSC run.
- **Outputs are written back next to the source h5 files.** With `--pipeline-h5`,
  the corrected matrix is saved as
  `h5_files/<SampleID>/sample_cellbender_filtered.h5` — *alongside* the Cell Ranger
  files, nothing is overwritten. (The full working outputs — report, checkpoint,
  posterior, optional `.h5ad` — go in a separate `cellbender_samples/` folder.)
- **The R pipeline then selects whether to load it.** Script 01 has a toggle:
  set `USE_CELLBENDER <- TRUE` and it ingests `sample_cellbender_filtered.h5`
  instead of `sample_filtered_feature_bc_matrix.h5`. Set `RUN_DECONTX <- FALSE`
  so ambient RNA isn't corrected twice.

So the flow is: **run CellBender in the terminal → files land in `h5_files/` →
flip one switch in Script 01.**

## 1. Environment (one time)

CellBender needs its own env (`environment_cellbender.yml` here); it also includes
scanpy/anndata for the optional `.h5ad` conversion, so this one env runs the whole
script.

```bash
conda env create -f environment_cellbender.yml
conda activate cellbender_env
cellbender remove-background --help     # sanity check
```

This installs the **stable release v0.4.0** (same outputs as v0.3.2, Python
3.10–3.12). For a **GPU** build of PyTorch, install the CUDA wheel matching your
driver *instead of* the plain `torch` line — see the comment in the yml, e.g.
`pip install torch --index-url https://download.pytorch.org/whl/cu121`.
Do **not** install `@master` (dev) or v0.3.1 (redacted for a bug).

## 2. Inputs

- A metadata `.xlsx` with at least a **`SampleID`** column (matching the folder
  names under `h5_files/`). Any extra columns are carried onto the `.h5ad` `.obs`.
- Optional per-sample columns **`expected_cells`** and **`total_droplets`** — if
  present they override the global defaults for that sample (blank = default).
- Per-sample raw matrices staged by `copy_sc_data.sh`:
  `h5_files/<SampleID>/sample_raw_feature_bc_matrix.h5`.

## 3. Run it (external terminal)

From the project's `r_process` directory (where `h5_files/` lives):

```bash
conda activate cellbender_env

CB_METADATA=Nr4a1_s17_metadata.xlsx \
CB_H5_DIR=h5_files \
python /path/to/cellbender/process_samples_raw_cellbender.py --cuda --pipeline-h5
```

- `--cuda` uses the GPU (auto-detected if omitted); `--cpu` forces CPU (slow).
- `--pipeline-h5` writes `h5_files/<SampleID>/sample_cellbender_filtered.h5` for
  Script 01. Omit it if you only want the scanpy-side `.h5ad`/report outputs.
- `--dry-run` prints the per-sample commands without running anything.
- Already-finished samples are skipped, so you can re-run to resume.

Everything is overridable by env var or flag (run with `--help`):
paths, `--sample-col`, `--raw-h5`, `--expected-cells`, `--total-droplets`,
`--epochs`, `--fpr`, `--learning-rate`, `--extra-args`, `--cellbender <path>`.

**Outputs**
- `cellbender_samples/<sample>/` — `*_cellbender.h5`, `*_cellbender_filtered.h5`,
  `*_cellbender.h5ad` (if scanpy present), `*.html` report, `*_ckpt.tar.gz`.
- `h5_files/<SampleID>/sample_cellbender_filtered.h5` (only with `--pipeline-h5`).
- `cellbender_run_log.txt` — per-sample run log + success/fail tally.

**Always inspect the per-sample HTML report** before trusting the output
(cell-calling, learning curve, removed counts).

## 4. Use it in the R pipeline (Script 01)

In `01_process_data.R`:

```r
USE_CELLBENDER <- TRUE      # ingest the CellBender matrices
CELLBENDER_H5  <- "sample_cellbender_filtered.h5"   # (default; leave as is)
RUN_DECONTX    <- FALSE     # CellBender already removed ambient RNA
```

Script 01 then reads `h5_files/<SampleID>/sample_cellbender_filtered.h5` per
sample; everything downstream (QC, doublets, merge, annotation) is unchanged.
To go back to the Cell Ranger + DecontX path, set `USE_CELLBENDER <- FALSE`
(and `RUN_DECONTX <- TRUE`).

## Files here

| file | purpose |
|------|---------|
| `process_samples_raw_cellbender.py` | portable per-sample CellBender runner |
| `environment_cellbender.yml`        | conda env (CellBender v0.4.0 + PyTorch + scanpy) — runs the whole script |
| `metadata_template.xlsx`             | example metadata layout (needs a `SampleID` column; add optional `expected_cells`/`total_droplets` for per-sample params) |
