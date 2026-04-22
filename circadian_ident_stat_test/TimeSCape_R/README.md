# TimeSCape v0.2 — R  (Seurat + SingleCellExperiment)

Circadian rhythm detection pipeline for single-cell RNA-seq data.
Direct translation of the MATLAB v0.2 pipeline, accepting **Seurat** (v4/v5) and **SingleCellExperiment** (Bioconductor) objects as input.
Includes a Shiny web GUI as the equivalent of the MATLAB `TimeSCape_GUI.m`, plus a full pathway-level circadian pipeline (GOBP gene sets → AUCell scoring → pathway cosinor → GRN time series).

---

## Installation

### 1 — Core CRAN packages

Install all at once. Use `dependencies = TRUE` to avoid missing sub-dependencies:

```r
install.packages(c(
  "Seurat",           # Seurat input format (v4 or v5)
  "minpack.lm",       # Levenberg-Marquardt NLS solver — REQUIRED
  "future",           # parallel workers
  "future.apply",     # future_lapply (block parallelism)
  "pheatmap",         # heatmap rendering
  "ggplot2",          # gene / pathway plots
  "gridExtra",        # multi-panel plot grids
  "openxlsx",         # Excel output for pathway results
  "igraph",           # GRN network construction
  "ggraph",           # GRN network visualisation
  "cowplot",          # GRN panel assembly
  "msigdbr",          # MSigDB gene set download (GOBP)
  "shiny",            # GUI (optional, only for app.R)
  "bslib",            # GUI theme
  "rhandsontable"     # editable ZT table in GUI
), dependencies = TRUE)
```

### 2 — Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "AUCell",                  # pathway scoring (required for pathway pipeline)
  "BiocParallel",            # AUCell parallelism
  "SingleCellExperiment",    # SCE input format (optional)
  "SummarizedExperiment"     # SCE input format (optional)
))
```

> **Verify installation before running:**
>
> ```r
> required_cran <- c("minpack.lm", "future", "future.apply", "pheatmap",
>                    "ggplot2", "gridExtra", "openxlsx", "igraph", "ggraph",
>                    "cowplot", "msigdbr", "shiny", "bslib", "rhandsontable")
> required_bioc <- c("AUCell", "BiocParallel")
> missing_cran  <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
> missing_bioc  <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
> if (length(missing_cran) > 0) install.packages(missing_cran, dependencies = TRUE)
> if (length(missing_bioc) > 0) BiocManager::install(missing_bioc)
> message("All packages ready.")
> ```

### 3 — Source the functions

```r
source("R/estimate_phaseR.R")
source("R/run_timescape.R")
source("R/generate_heatmap.R")
source("R/plot_gene.R")
source("R/pathway_circadian.R")   # pathway pipeline (GOBP + AUCell + GRN)
```

Or simply open `app.R` in RStudio and click **Run App** — it sources everything automatically.

---

## Quick Start — Command Line

```r
library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)

set.seed(123)

# Source all pipeline functions
source("R/estimate_phaseR.R")
source("R/run_timescape.R")
source("R/generate_heatmap.R")
source("R/plot_gene.R")
source("R/pathway_circadian.R")

# Load your Seurat object
obj <- readRDS("my_data.rds")

# Auto-parse ZT times (expects "ZT00", "ZT06", ... in the ZT column)
tmeta <- build_tmeta(obj, zt_col = "ZT_time_str")
print(tmeta)   # verify ZT_times; edit manually if your naming differs

# Run full gene-level circadian pipeline
results <- run_timescape(
  obj          = obj,
  celltype_col = "CellType",
  zt_col       = "ZT_time_str",
  tmeta        = tmeta,
  norm_str     = "logcounts",   # see Normalization section below
  n_workers    = 2L,            # parallel workers (2 = good balance)
  outdir       = "TimeSCape_output"
)
```

---

## Quick Start — Shiny GUI

```r
shiny::runApp("app.R")
```

Or in RStudio: open `app.R` → click **Run App** (top-right of the editor pane).

Opens in your default browser. Workflow:

| Step | Action |
|------|--------|
| ① | Browse and load a `.rds` Seurat file |
| ② | Pick cell-type column, ZT column, optional group column |
| ② | Click **Build / Preview ZT Table** → verify numeric ZT hours |
| ③ | Choose normalization, period, workers; click **▶ Run Analysis** or **▶▶ All Types** |
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
| `n_workers` | `1L` | Parallel workers. `1` = sequential (safest); `2` = 2.5× faster for sparse data |
| `rm_low_conf` | `TRUE` | Write confident-only CSV subsets |
| `plot_heat` | `TRUE` | Save heatmap PNG after each cell type |
| `custom_genelist` | `NULL` | Restrict to these genes only |
| `custom_celltype` | `NULL` | Restrict to these cell types only |
| `group_col` | `NULL` | **Optional second grouping column** (see below) |
| `custom_group` | `NULL` | Filter to specific group values |
| `custom_zt` | `NULL` | Filter to specific ZT time points (e.g. `c("ZT00","ZT12")`) |
| `outdir` | `getwd()` | Root output directory |

---

## Normalization (`norm_str`)

| Value | When to use |
|-------|-------------|
| `"lib_size"` | **Recommended for standard scRNA-seq.** Per-cell library-size normalisation to 10 000 counts + `log1p`. Identical to MATLAB `pkg.norm_libsize(X, 1e4)`. Requires a non-empty `counts` slot. |
| `"logcounts"` | **Required for decontX / SCTransform outputs.** Use when the `counts` slot is empty but the `data` slot holds pre-normalised expression (decontX stores cleaned data here). Also correct after Seurat `NormalizeData()` or SCE `logNormCounts()`. |
| `"none"` | Use raw counts as-is. Only when data is already normalised externally and stored in `counts`. |

> **decontX users**: use `norm_str = "logcounts"`. decontX stores cleaned, normalised expression in the Seurat `data` slot; the `counts` slot is empty. Using `"lib_size"` would silently return 0 for all genes.

---

## Parallel Workers (`n_workers`)

`n_workers = 1` (sequential) is the safest option — no data copies, minimum RAM.

`n_workers = 2` gives approximately 2.5× speedup for sparse scRNA-seq matrices. Workers receive only the per-cell-type sparse slice (not the full Seurat object), so RAM overhead is modest. Requires `future` and `future.apply` packages.

```r
# Two workers — validated at ~2.5x speedup, ~2x RAM vs sequential
results <- run_timescape(obj, ..., n_workers = 2L)
```

> **Dense matrices (decontX, MAGIC-imputed):** prefer `n_workers = 1`. Dense matrices are large even when sliced per cell type, and serialising them to two workers doubles peak RAM.

---

## Second Grouping Variable (`group_col`)

When all cancer stages (or replicates, treatments, …) are in a single object, use `group_col` to run the analysis independently for every **(cell type × group)** combination without splitting the object:

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

Output structure mirrors cell-type dirs:
```
TimeSCape_output/
  CD8_T_cells_Early/
  CD8_T_cells_Late/
  all_cell_types_period_24_summary_results.csv
```

---

## Pathway Circadian Pipeline (`pathway_circadian.R`)

After gene-level analysis, run the full pathway pipeline:

```r
# 1. Download deduplicated GOBP gene sets from MSigDB
genesets <- pull_genesets(
  collection  = "C5",
  subcategory = "GO:BP",
  species     = "Mus musculus",   # or "Homo sapiens"
  min_size    = 10L,
  max_size    = 500L,
  deduplicate = TRUE   # removes POSITIVE/NEGATIVE/REGULATION_OF_X redundancy
)

# 2. Score cells with AUCell (cache result — takes several minutes for >100k cells)
auc_cache <- "TimeSCape_output/auc_matrix.rds"
if (file.exists(auc_cache)) {
  auc_matrix <- readRDS(auc_cache)
} else {
  auc_matrix <- auc_score_cells(obj, genesets, use_norm = TRUE, auc_max_rank = 0.05)
  saveRDS(auc_matrix, auc_cache)
}

# 3. Pathway-level cosinor (same model as gene-level, but per-cell AUC scores)
path_results <- pathway_cosinor(
  auc_mat      = auc_matrix,
  meta         = obj@meta.data,
  celltype_col = "CellType",
  zt_col       = "ZT_time_str",
  tmeta        = tmeta,
  target_ct    = "CD8+ T cells"
)

# 4. Write Excel results (two sheets: All + Confident)
write_pathway_results(path_results, "TimeSCape_output/CD8_T_cells_pathway_circadian.xlsx",
                      celltype = "CD8+ T cells")

# 5. Plot top pathways
plot_pathway_single(auc_matrix, path_results, obj@meta.data,
                    celltype_col = "CellType", zt_col = "ZT_time_str",
                    tmeta = tmeta, target_ct = "CD8+ T cells",
                    target_pathway = "GOBP_T_CELL_ACTIVATION")

# 6. GRN time series — correlation network per ZT, coloured by gene type
plot_grn_timeseries(
  obj           = obj,
  circ_genes    = conf_genes,
  pathway_genes = genesets[["GOBP_T_CELL_ACTIVATION"]],
  meta          = obj@meta.data,
  celltype_col  = "CellType",
  zt_col        = "ZT_time_str",
  tmeta         = tmeta,
  target_ct     = "CD8+ T cells",
  cor_thresh    = 0.20,
  p_thresh      = 0.05
)
```

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
| `*_pathway_circadian.xlsx` | Pathway cosinor results (All + Confident sheets) |
| `*_GRN_timeseries_*.png` | GRN network panels, one per ZT time point |

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
| `pvalue_adj` | BH-corrected F-test p-value (denominator = all genes including zero-expression) |
| `Sine_corr` | Pearson r (per-ZT means vs fitted cosine) |
| `pvalue_corr` | Pearson correlation p-value (raw) |
| `pvalue_adj_corr` | BH-corrected correlation p-value |

Rows sorted: `pvalue_adj_corr ↑ → pvalue_adj ↑ → Acrophase_24 ↑ → Abs_Amp ↓`

**Confident criterion**: `pvalue < 0.05` **AND** `pvalue_corr < 0.05` (raw, not BH-adjusted).

---

## ZT Metadata — `build_tmeta()`

Parses numeric ZT hours from a string metadata column. Supports:
`"ZT00"`, `"ZT06"`, `"zt12"`, `"ZT_06"`, `"ZT 3"`, `"0"`, `"6"` …

```r
tmeta <- build_tmeta(obj, zt_col = "ZT_time_str")
# Returns data.frame: zt_str, ZT_times
# Edit manually if auto-parsing fails:
tmeta$ZT_times[tmeta$zt_str == "custom_label"] <- 6
```

`build_tmeta_from_seurat()` is kept as a backward-compatible alias.

---

## Individual Gene Plot

```r
p <- plot_gene_single(
  tmeta        = tmeta,
  cust_cells   = "CD8_T_cells",   # sanitised cell-type name (matches output dir)
  cust_gene    = "Nr1d1",
  period12     = FALSE,
  print_scdata = TRUE,
  sce          = obj,
  celltype_col = "CellType",
  zt_col       = "ZT_time_str",
  use_violin   = FALSE,   # FALSE = jittered dots; TRUE = violin density
  outdir       = "TimeSCape_output/CD8_T_cells"
)
print(p)
```

---

## Batch Gene Plots

```r
save_batch_plots(
  tmeta      = tmeta,
  cust_cells = "CD8_T_cells",
  plot_type  = 1,    # 1 = confident, 2 = non-confident, 3 = classical circadian
  period12   = FALSE,
  outdir     = "TimeSCape_output/CD8_T_cells"
)
```

---

## Troubleshooting

**"0 genes tested" / all genes return NA**

Run `diagnose_phaseR()` on a few genes to reveal the underlying error:
```r
diagnose_phaseR(obj, celltype_col = "CellType", zt_col = "ZT_time_str",
                tmeta = tmeta, target_ct = "CD8+ T cells",
                norm_str = "logcounts", n = 5L)
```
Common causes:
- `minpack.lm` not installed → `install.packages("minpack.lm")`
- `norm_str = "lib_size"` with empty `counts` slot → switch to `norm_str = "logcounts"`
- Fewer than 4 ZT time points for that cell type

**pvalue_adj values differ slightly from MATLAB**

This is expected if the gene count differs between runs. From v0.2 onward, BH correction includes all genes (zero-expression genes treated as p = 1), matching MATLAB's denominator. Raw p-values are unaffected.

**Parallel workers: S4 dispatch error (`[` not subsettable)**

Add `future.packages = c("Matrix", "minpack.lm")` (already handled internally). If you still see this, update to the current `run_timescape.R`.

---

## Differences from MATLAB v0.2

| Aspect | MATLAB v0.2 | R v0.2 |
|--------|-------------|--------|
| Input format | `SingleCellExperiment` (scGEAtoolbox) | `Seurat` (v4/v5) **or** `SingleCellExperiment` — auto-detected |
| NLS solver | Trust-Region (`fit()`) | Levenberg-Marquardt (`minpack.lm::nlsLM`) |
| Parallelism | `parfor` (implicit pool) | `future.apply::future_lapply()` + `plan()` |
| Heatmap | `imagesc` + custom colormap | `pheatmap` |
| Gene plots | MATLAB `plot` / `violinplot` | `ggplot2` |
| GUI | `uicontrol` figures | Shiny web app (`app.R`) |
| Group column | Not supported | `group_col` parameter |
| Pathway pipeline | Not included | `pathway_circadian.R` (GOBP + AUCell + cosinor + GRN) |
| BH denominator | All genes (incl. zero-expr.) | All genes (incl. zero-expr. as p=1) — matched from v0.2 |
| Numerical agreement | — | < 0.1 % difference vs MATLAB (validated 96.8% gene overlap) |

---

## File Structure

```
TimeSCape_R/
├── app.R                        # Shiny GUI — shiny::runApp("app.R") or RStudio Run App
├── README.md                    # This file
├── INSTALL.md                   # Extended installation notes
├── DESCRIPTION                  # R package metadata
└── R/
    ├── estimate_phaseR.R        # Cosinor NLS fitting + F-test + Pearson corr
    ├── run_timescape.R          # Main pipeline (cell-type × group loop, parallel)
    ├── generate_heatmap.R       # pheatmap z-score heatmap
    ├── plot_gene.R              # plot_gene_single() + save_batch_plots()
    └── pathway_circadian.R      # GOBP gene sets + AUCell + pathway cosinor + GRN
```

---

## Version

**v0.2** — April 2026
