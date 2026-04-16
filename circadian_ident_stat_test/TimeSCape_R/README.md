# TimeSCape v0.2 — R  (Seurat + SingleCellExperiment)

Circadian rhythm detection pipeline for single-cell RNA-seq data.
Direct translation of the MATLAB v0.2 pipeline, accepting **Seurat** (v4/v5) and **SingleCellExperiment** (Bioconductor) objects as input.
Includes a Shiny web GUI as the equivalent of the MATLAB `TimeSCape_GUI.m`.

---

## Installation

### 1 — R packages

Install all CRAN packages at once. Use `dependencies = TRUE` to avoid missing sub-dependencies (known issue with `minpack.lm` and `rhandsontable` without it):

```r
install.packages(c(
  "Seurat",           # Seurat input format (v4 or v5)
  "minpack.lm",       # Levenberg-Marquardt NLS solver
  "future",           # parallel workers
  "future.apply",     # future_lapply (block parallelism)
  "pheatmap",         # heatmap rendering
  "ggplot2",          # gene plots
  "shiny",            # GUI (optional, only for app.R)
  "bslib",            # GUI theme
  "rhandsontable"     # editable ZT table in GUI
), dependencies = TRUE)
```

For **SingleCellExperiment** input (Bioconductor), also install:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment"))
```

> **Tip — verify installation:** After installing, run the block below to catch any missing packages before your analysis:
>
> ```r
> required <- c("minpack.lm", "future", "future.apply", "pheatmap",
>               "ggplot2", "shiny", "bslib", "rhandsontable")
> missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
> if (length(missing) > 0) {
>   install.packages(missing, dependencies = TRUE)
> } else {
>   message("All packages installed.")
> }
> ```

### 2 — Source the functions

```r
source("R/estimate_phaseR.R")
source("R/run_timescape.R")
source("R/generate_heatmap.R")
source("R/plot_gene.R")
```

Or simply open `app.R` in RStudio and click **Run App** — it sources everything automatically.

---

## Quick Start — Command Line

```r
library(Seurat)   # or library(SingleCellExperiment)
library(future)

# Optional: parallel workers (recommended for > 5 000 genes)
plan(multisession, workers = 4)

# Load your object — Seurat or SingleCellExperiment both work
obj <- readRDS("my_data.rds")

# Preview metadata columns and ZT strings
head(obj@meta.data)           # Seurat
# or: head(as.data.frame(colData(obj)))   # SingleCellExperiment

# Auto-parse ZT times (expects "ZT00", "ZT06", ... in the ZT column)
tmeta <- build_tmeta(obj, zt_col = "ZT_time_str")
print(tmeta)   # verify ZT_times; edit manually if your naming differs

# Run full pipeline
results <- run_timescape(
  obj          = obj,
  celltype_col = "CellType",      # column with cell-type labels
  zt_col       = "ZT_time_str",   # column with ZT time strings
  tmeta        = tmeta,
  norm_str     = "lib_size",      # recommended
  outdir       = "TimeSCape_output"
)
```

---

## Quick Start — Shiny GUI

```r
shiny::runApp("app.R")
```

Opens in your browser. Workflow mirrors the MATLAB GUI:

| Step | Action |
|------|--------|
| ① | Browse and load a `.rds` Seurat file |
| ② | Pick cell-type column, ZT column, optional group column |
| ② | Click **Build / Preview ZT Table** → verify numeric ZT hours |
| ③ | Choose normalization, period; click **▶ Run Analysis** or **▶▶ All Types** |
| ④ | Click **Generate Heatmap** |
| ⑤ | Select gene → **Plot Single Gene**; or **Save (Batch)** for all genes |
| — | Download PNG / PDF / SVG from the plot header |

---

## Key Parameters — `run_timescape()`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `obj` | — | Seurat (v4/v5) **or** SingleCellExperiment object — auto-detected |
| `celltype_col` | `"cell_type"` | Metadata column with cell-type labels |
| `zt_col` | `"ZT_str"` | Metadata column with ZT strings (`"ZT00"`, `"ZT06"`, …) |
| `tmeta` | `NULL` | data.frame from `build_tmeta()`. Auto-built if NULL. |
| `norm_str` | `"lib_size"` | See **Normalization** section below |
| `period12` | `FALSE` | `TRUE` = 12-hr period; `FALSE` = 24-hr (default) |
| `rm_low_conf` | `TRUE` | Write confident-only CSV subsets |
| `plot_heat` | `TRUE` | Save heatmap PNG after each cell type |
| `custom_genelist` | `NULL` | Restrict to these genes only |
| `custom_celltype` | `NULL` | Restrict to these cell types only |
| `group_col` | `NULL` | **Optional second grouping column** (see below) |
| `custom_group` | `NULL` | Filter to specific group values |
| `outdir` | `getwd()` | Root output directory |

---

## Normalization (`norm_str`)

| Value | Behaviour |
|-------|-----------|
| `"lib_size"` | Per-cell library-size normalisation to 10 000 counts + `log1p`. Identical to MATLAB `pkg.norm_libsize(X, 1e4)`. **Recommended.** Safe across cancer stages and replicates because each cell's factor depends only on its own total count sum. |
| `"logcounts"` | Use a pre-computed normalised slot: Seurat `NormalizedData` (`data` slot) or SCE `logcounts` assay. Use if you already ran `NormalizeData()` / `SCTransform` (Seurat) or `scuttle::logNormCounts()` (SCE). |
| `"none"` | Use raw counts as-is from the `counts` slot. Use when data is already normalised externally and stored in the counts slot, or when you explicitly want no transformation. |

Normalisation happens **inside the cell-type loop** (or cell-type × group loop) so only one dense submatrix is in RAM at a time — memory-efficient even for large datasets.

---

## Second Grouping Variable (`group_col`)

When all cancer stages (or replicates, treatments, …) are in a single Seurat object, use `group_col` to run the analysis independently for every **(cell type × group)** combination without splitting the object manually.

```r
results <- run_timescape(
  obj          = obj,
  celltype_col = "CellType",
  zt_col       = "ZT_time_str",
  group_col    = "tumor_stage",    # e.g. "Early", "Intermediate", "Late"
  norm_str     = "lib_size",
  outdir       = "TimeSCape_output"
)
```

**Output structure:**

```
TimeSCape_output/
  Hepatocytes_early/
    Hepatocytes_early_period_24_circadian_analysis_all.csv
    Hepatocytes_early_period_24_circadian_analysis_confident.csv
    Hepatocytes_early_period_24_heatmap_strict.png
    ...
  Hepatocytes_mid/
    ...
  Hepatocytes_late/
    ...
  T_cells_early/
    ...
  all_cell_types_period_24_summary_results.csv   ← Group column included
```

The returned list is keyed by the combo name:

```r
results[["Hepatocytes_early"]]$T1   # stats table
results[["Hepatocytes_early"]]$T2   # per-ZT means
```

In the **Shiny GUI**, select the group column in step ② and a group value picker appears in the Gene Explorer — no code changes needed.

---

## Output Files (per cell type / combo)

| File | Contents |
|------|----------|
| `*_circadian_analysis_all.csv` | All genes tested — stats table (see columns below) |
| `*_circadian_analysis_confident.csv` | Genes with F-test **and** corr p < 0.05 (raw) |
| `*_circadian_ZTs_mean.csv` | Per-ZT mean expression (rows = genes) |
| `*_circadian_ZTs_mean_normalized.csv` | ZT00-normalised means |
| `*_circadian_ZTs_mean_confident.csv` | Per-ZT means, confident genes only |
| `*_circadian_ZTs_mean_normalized_confident.csv` | Normalised, confident only |
| `*_heatmap_strict.png` | Blue→white→red z-score heatmap, confident genes |
| `*_heatmap_genes.csv` | Gene table used to draw the heatmap |
| `all_cell_types_*_summary_results.csv` | Cross-cell-type / cross-group summary |

### Stats table columns

| Column | Definition |
|--------|-----------|
| `Genes` | Gene name |
| `Amp` | Cosine amplitude (signed; negative = peak shifted by 12 h) |
| `Abs_Amp` | \|Amplitude\| |
| `Mesor` | Midline estimating statistic of rhythm (mean expression level) |
| `Acrophase` | Raw peak time from NLS fit (hrs) |
| `Acrophase_24` | Peak time wrapped to [0, 24) |
| `Period` | 12 or 24 hrs |
| `pvalue` | F-test p-value |
| `pvalue_adj` | BH-corrected F-test p-value |
| `Sine_corr` | Pearson r (per-ZT means vs fitted cosine) |
| `pvalue_corr` | Pearson correlation p-value (raw) |
| `pvalue_adj_corr` | BH-corrected correlation p-value |

Rows sorted: `pvalue_adj_corr ↑ → pvalue_adj ↑ → Acrophase_24 ↑ → Abs_Amp ↓`

**Confident criterion**: `pvalue < 0.05` **AND** `pvalue_corr < 0.05` (raw, not BH-adjusted).

---

## ZT Metadata — `build_tmeta()`

Parses numeric ZT hours from a string metadata column. Works with both Seurat and SingleCellExperiment. Supports:
`"ZT00"`, `"ZT06"`, `"zt12"`, `"ZT_06"`, `"ZT 3"`, `"0"`, `"6"` …

```r
tmeta <- build_tmeta(obj, zt_col = "ZT_time_str")
# Returns data.frame with columns: zt_str, ZT_times
# Edit ZT_times manually if auto-parsing fails for your naming convention:
tmeta$ZT_times[tmeta$zt_str == "custom_label"] <- 6
```

`build_tmeta_from_seurat()` is kept as a backward-compatible alias.

---

## Individual Gene Plot

```r
p <- plot_gene_single(
  tmeta        = tmeta,
  cust_cells   = "Hepatocytes",   # must match directory name (sanitised)
  cust_gene    = "Nr1d1",
  period12     = FALSE,
  print_scdata = TRUE,            # overlay single-cell data
  sce          = seu,             # Seurat object for overlay
  celltype_col = "cell_type",
  zt_col       = "ZT_str",
  use_violin   = TRUE,            # TRUE = violin, FALSE = dots scatter
  outdir       = "TimeSCape_output/Hepatocytes"
)
print(p)   # ggplot2 object — fully customisable
```

When `group_col` is active, pass the combo directory:

```r
p <- plot_gene_single(
  ...,
  cust_cells = "Hepatocytes_early",
  outdir     = "TimeSCape_output/Hepatocytes_early"
)
```

---

## Batch Gene Plots

```r
save_batch_plots(
  tmeta      = tmeta,
  cust_cells = "Hepatocytes",
  plot_type  = 1,   # 1 = confident, 2 = non-confident, 3 = classical circadian
  period12   = FALSE,
  outdir     = "TimeSCape_output/Hepatocytes"
)
# Saves one PNG per gene to outdir/Hepatocytes_period_24_plots_confident/
```

---

## Differences from MATLAB v0.2

| Aspect | MATLAB v0.2 | R v0.2 |
|--------|-------------|--------|
| Input format | `SingleCellExperiment` (scGEAtoolbox) | `Seurat` (v4/v5) **or** `SingleCellExperiment` — auto-detected |
| NLS solver | Trust-Region (`fit()`) | Levenberg-Marquardt (`minpack.lm::nlsLM`) |
| Parallelism | `parfor` (implicit pool) | `future.apply::future_lapply()` + `plan()` |
| Heatmap | `imagesc` + custom colormap | `pheatmap` |
| Gene plots | MATLAB `plot` / `violinplot` | `ggplot2` |
| GUI | `uicontrol` figures | Shiny web app |
| Group column | Not supported | `group_col` parameter |
| Numerical agreement | — | < 0.1 % difference vs MATLAB |

---

## File Structure

```
TimeSCape_R/
├── app.R                    # Shiny GUI — run with shiny::runApp("app.R")
├── README.md                # This file
├── INSTALL.md               # Extended installation notes
├── DESCRIPTION              # R package metadata
└── R/
    ├── estimate_phaseR.R    # Cosinor NLS fitting + F-test + Pearson corr
    ├── run_timescape.R      # Main pipeline (cell-type × group loop)
    ├── generate_heatmap.R   # pheatmap z-score heatmap
    └── plot_gene.R          # plot_gene_single() + save_batch_plots()
```

---

## Version

**v0.2** — April 2026
