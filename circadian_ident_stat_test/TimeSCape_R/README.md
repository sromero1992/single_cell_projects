# TimeSCape v0.2 — R  (Seurat + SingleCellExperiment)

Circadian rhythm detection pipeline for single-cell RNA-seq data.
Direct translation of the MATLAB v0.2 pipeline, extended with a full
pathway-level circadian screen (two complementary approaches).
Accepts **Seurat** (v4/v5) and **SingleCellExperiment** objects.

---

## Repository layout

```
TimeSCape_R/
├── README.md                          # This file
├── INSTALL.md                         # Extended installation notes
├── app.R                              # Shiny GUI (optional)
│
├── R/                                 # Source functions — do not edit for routine use
│   ├── estimate_phaseR.R              # Cosinor NLS fitting + F-test + Pearson corr
│   ├── run_timescape.R                # Gene-level pipeline (cell-type loop, parallel)
│   ├── generate_heatmap.R             # pheatmap z-score heatmap
│   ├── plot_gene.R                    # Single-gene plot (returns ggplot object)
│   └── pathway_circadian.R            # AUCell scoring + pathway cosinor + plots
│
├── run_timescape_test.R               # ← Single cell-type debugger (CD8+ T cells)
│                                      #   Use this to test parameters before a full run
│
├── run_full_pipeline_early.R          # ← Production: all cell types, early stage
├── run_full_pipeline_intermediate.R   #   all cell types, intermediate stage
└── run_full_pipeline_advanced.R       #   all cell types, advanced stage
```

---

## Three-step pipeline

Each production script runs three sequential analyses **per cell type**:

### Step A — Gene-level cosinor

`run_timescape()` fits a 24-hr (or 12-hr) cosine model to every gene using
NLS (Levenberg-Marquardt) + F-test + Pearson correlation.
Genes passing both p < 0.05 thresholds are **confident circadian genes**.

Outputs: all CSVs + confident-genes Excel + z-score heatmap PNG.
*(No individual gene plots — use `run_timescape_test.R` for those.)*

### Step B — Approach A: enricher() ORA → pathway cosinor

`clusterProfiler::enricher()` tests which pathways are significantly
enriched among the confident circadian genes (BH-corrected ORA).
Cells are then scored for each enriched pathway (AUCell or AddModuleScore),
and `pathway_cosinor()` identifies which of those pathway scores oscillate.

This answers: *"Which known pathways contain circadian-regulated genes?"*

### Step C — Approach B: all-pathway circadian screen

**All** size-filtered pathways in the selected databases (KEGG, Reactome,
GO:BP, GO:MF) are scored on the subsetted cells, and `pathway_cosinor()`
identifies oscillating gene programs — without any gene pre-filter.

This answers: *"Which gene programs oscillate rhythmically, even if the
individual member genes are not significantly circadian on their own?"*

Approaches A and B are complementary:
A finds pathways *enriched in* individually-circadian genes;
B finds pathways whose *collective expression* is rhythmic.

---

## Installation

### 1 — Core CRAN packages

```r
install.packages(c(
  "Seurat",           # Seurat input format (v4 or v5)
  "minpack.lm",       # Levenberg-Marquardt NLS solver — REQUIRED
  "future",           # parallel workers
  "future.apply",     # future_lapply (block parallelism)
  "pheatmap",         # heatmap rendering
  "ggplot2",          # gene / pathway plots
  "gridExtra",        # multi-panel plot grids
  "openxlsx",         # Excel output
  "msigdbr",          # MSigDB gene set download
  "dplyr",            # data wrangling
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
  "clusterProfiler",       # enricher() ORA (Step B)
  "AUCell",                # rank-based pathway scoring
  "BiocParallel",          # AUCell parallelism
  "SingleCellExperiment",  # SCE input format (optional)
  "SummarizedExperiment"   # SCE input format (optional)
))
```

### 3 — Verify and source

```r
# Quick dependency check
required <- c("minpack.lm","future","future.apply","pheatmap","ggplot2",
              "gridExtra","openxlsx","msigdbr","dplyr","clusterProfiler","AUCell")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) message("Missing: ", paste(missing, collapse=", "))

# Source all functions
src <- "C:/path/to/TimeSCape_R/R"
source(file.path(src, "estimate_phaseR.R"))
source(file.path(src, "run_timescape.R"))
source(file.path(src, "generate_heatmap.R"))
source(file.path(src, "plot_gene.R"))
source(file.path(src, "pathway_circadian.R"))
```

---

## Quick start — production run

Open one of the three stage scripts and set the four paths at the top:

```r
# run_full_pipeline_early.R  (lines 39-43)
stage    <- "early"
base_dir <- "Z:/your/data/directory"
src_path <- "C:/path/to/TimeSCape_R/R"
rds_file <- "your_seurat_object_early.rds"
out_dir  <- file.path(base_dir, sprintf("TimeSCape_output_%s", stage))
```

Then source or `Rscript` the file. It will:
1. Load the RDS
2. Merge cell-type subtypes (edit the merging block in Section 4 as needed)
3. Build TERM2GENE from msigdbr once
4. Loop over every cell type → Step A → Step B → Step C
5. Write a `run_log_{stage}.csv` summary

---

## Quick start — single cell type (debug / exploration)

Use `run_timescape_test.R` to iterate on one cell type before committing
to a full run. It includes the same three-step pipeline plus:
- Interactive top-6 gene grid
- Single-gene explorer block (Section 6d) — tweak `explore_gene` and re-run

---

## Key parameters (production scripts)

### Scoring method

```r
use_aucell <- TRUE    # AUCell: rank-based, robust (recommended)
                      # FALSE = AddModuleScore: faster, ctrl=5 for prototyping
```

AUCell ranks genes per cell once and evaluates each pathway as a fast set
lookup — much faster than AddModuleScore when scoring thousands of pathways.
Switch to `FALSE` only if AUCell is not installed.

### Database selection

```r
use_kegg     <- TRUE    # ~200 curated metabolic/signalling pathways
use_reactome <- TRUE    # ~1,500 fine-grained pathway hierarchy
use_gobp     <- TRUE    # ~7,500 biological process terms (large; set FALSE for speed)
use_gomf     <- FALSE   # ~1,000 molecular function terms (often redundant with GO:BP)
```

### Gene set size filter

```r
min_gs_size <- 10L    # pathways with fewer measured genes are skipped
max_gs_size <- 500L   # very large gene sets are skipped (genome-wide noise)
```

Applied after intersecting with genes actually measured in your object.
Typically reduces ~9,000 raw pathways to ~5,000 scoreable gene sets.

### Enricher ORA (Step B)

```r
enrich_pval   <- 0.05   # nominal p-value cutoff
enrich_padj   <- 0.20   # BH-adjusted cutoff (relaxed; pathway cosinor is the final filter)
enrich_min_gs <- 10L
enrich_max_gs <- 500L
```

### Cosinor thresholds (Steps B & C)

```r
cosinor_pval      <- 0.05   # F-test p-value
cosinor_pval_corr <- 0.05   # Pearson correlation p-value
```

Both must pass for a pathway to be called **oscillating**.

### Cell subtype merging

Section 4 of each production script collapses sub-annotations into clean
groups before the loop. Edit to match your annotation scheme:

```r
data$cell_types_merged[data$sub_cell_types2 %in%
  c("CD8+ T cells", "Cyc. CD8+ T cells", "TCF7+ CD8+ T cells")] <- "CD8+ T cells"
```

---

## Output structure (per cell type per stage)

```
TimeSCape_output_early/
└── CD8_T_cells/
    │
    │  ── Step A: gene-level ──────────────────────────────────────────────
    ├── CD8_T_cells_period_24_circadian_analysis_all.csv
    ├── CD8_T_cells_period_24_circadian_analysis_confident.csv
    ├── CD8_T_cells_period_24_circadian_ZTs_mean.csv
    ├── CD8_T_cells_period_24_circadian_ZTs_mean_normalized.csv
    ├── CD8_T_cells_period_24_circadian_ZTs_mean_confident.csv
    ├── CD8_T_cells_period_24_circadian_ZTs_mean_normalized_confident.csv
    ├── CD8_T_cells_confident_genes.xlsx          ← clean sorted confident list
    ├── CD8_T_cells_period_24_heatmap_strict.png  ← z-score heatmap
    │
    │  ── Step B: Approach A ──────────────────────────────────────────────
    ├── CD8_T_cells_enricher_ORA.xlsx             ← enrichment results
    ├── CD8_T_cells_approachA_pathway_cosinor.xlsx
    │     Sheet "All_Pathways"        — cosinor stats + gene members
    │     Sheet "Oscillating_Pathways"— significant oscillators + gene members
    ├── CD8_T_cells_approachA_top_pathways.png    ← top-6 cosine grid
    ├── CD8_T_cells_enricherA_scores_aucell.rds   ← score cache (auto-reused)
    │
    │  ── Step C: Approach B ──────────────────────────────────────────────
    ├── CD8_T_cells_approachB_allpathways_cosinor.xlsx
    │     Sheet "All_Pathways"        — all scored pathways + gene members
    │     Sheet "Oscillating_Pathways"— significant oscillators + gene members
    ├── CD8_T_cells_approachB_top_pathways.png    ← top-6 cosine grid
    └── CD8_T_cells_allpathwaysB_scores_aucell.rds ← score cache (auto-reused)

run_log_early.csv   ← summary: cells, conf genes, enriched/oscillating counts per cell type
```

### Stats table columns (gene-level)

| Column | Definition |
|--------|------------|
| `Genes` | Gene symbol |
| `Amp` | Cosine amplitude (signed) |
| `Abs_Amp` | \|Amplitude\| |
| `Mesor` | Midline (mean expression) |
| `Acrophase` | Raw NLS peak time (hrs) |
| `Acrophase_24` | Peak time wrapped to [0, 24) |
| `Period` | 12 or 24 hrs |
| `pvalue` | F-test p-value |
| `pvalue_corr` | Pearson correlation p-value |
| `pvalue_corr_adj` | BH-corrected correlation p |

**Confident criterion**: `pvalue < 0.05` **AND** `pvalue_corr < 0.05`.

---

## Key parameters — `run_timescape()`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `obj` | — | Seurat (v4/v5) or SCE — auto-detected |
| `celltype_col` | `"cell_type"` | Metadata column with cell-type labels |
| `zt_col` | `"ZT_str"` | Metadata column with ZT strings (`"ZT00"`, `"ZT06"`, …) |
| `tmeta` | `NULL` | From `build_tmeta()`; auto-built if NULL |
| `norm_str` | `"lib_size"` | See **Normalization** below |
| `period12` | `FALSE` | `TRUE` = 12-hr; `FALSE` = 24-hr |
| `n_workers` | `1L` | Parallel workers (2 = ~2.5× faster for sparse data) |
| `rm_low_conf` | `TRUE` | Write confident-only CSV subsets |
| `plot_heat` | `TRUE` | Save heatmap after each cell type |
| `custom_celltype` | `NULL` | Restrict to specific cell types |
| `group_col` | `NULL` | Optional second grouping column (cancer stage, replicate…) |
| `outdir` | `getwd()` | Root output directory |

---

## Normalization (`norm_str`)

| Value | When to use |
|-------|-------------|
| `"lib_size"` | Standard scRNA-seq: per-cell normalisation to 10 000 + `log1p`. Requires non-empty `counts` slot. |
| `"logcounts"` | **decontX / SCTransform / `NormalizeData()` outputs.** Use when `counts` is empty and `data` holds pre-normalised values. |
| `"none"` | Data already normalised externally and stored in `counts`. |

> **decontX users**: always use `norm_str = "logcounts"`.

---

## ZT metadata — `build_tmeta()`

```r
tmeta <- build_tmeta(obj, zt_col = "ZT_time_str")
# Returns data.frame: zt_str, ZT_times  (sorted by hour)
# Edit manually if auto-parsing fails:
tmeta$ZT_times[tmeta$zt_str == "custom_label"] <- 6
```

Supports: `"ZT00"`, `"ZT06"`, `"zt12"`, `"ZT_06"`, `"ZT 3"`, `"0"`, `"6"` …

---

## Individual gene plot (interactive / explorer)

```r
p <- plot_gene_single(
  tmeta        = tmeta,
  cust_cells   = "CD8_T_cells",        # sanitised cell-type name
  cust_gene    = "Nr1d1",
  period12     = FALSE,
  print_scdata = TRUE,                  # overlay single-cell data
  sce          = obj,
  celltype_col = "cell_types_merged",
  zt_col       = "ZT_time_str",
  use_violin   = FALSE,                 # FALSE = jitter  |  TRUE = violin
  outdir       = "TimeSCape_output_early/CD8_T_cells"
)
print(p)
```

The `run_timescape_test.R` Section 6d provides a parameter block for
interactive exploration — set `explore_gene` and re-run the block.

---

## Shiny GUI

```r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio → **Run App**.

| Step | Action |
|------|--------|
| ① | Load `.rds` Seurat file or workspace variable |
| ② | Pick cell-type column, ZT column; define ZT times |
| ③ | Choose normalization, period, workers; click **▶ Run Analysis** or **▶▶ All Types** |
| ④ | Click **Generate Heatmap** (auto-generated for all types after Run All) |
| ⑤ | Select gene → **Plot Single Gene** |

---

## Troubleshooting

**"0 genes tested" / all genes return NA**
```r
diagnose_phaseR(obj, celltype_col = "cell_types_merged",
                zt_col = "ZT_time_str", tmeta = tmeta,
                target_ct = "CD8+ T cells", norm_str = "logcounts", n = 5L)
```
Common causes: `minpack.lm` not installed; `norm_str = "lib_size"` with empty
`counts` slot; fewer than 4 ZT time points for the cell type.

**Score cache is stale after changing databases or cell type merging**

Delete the `.rds` cache files in the cell type's output folder and re-run.
Cache filenames encode the method: `*_aucell.rds` or `*_ams.rds`.

**enricher() returns no pathways**

Relax `enrich_padj` (e.g. `0.50`) or check that `conf_genes` has at least
5 members. Very small cell types with few confident genes often fail the ORA.

**Parallel workers: S4 dispatch error**

Use `n_workers = 1L`. Dense matrices (decontX, MAGIC) double peak RAM when
serialised to workers; sequential mode is safer and only ~2× slower.

**Windows MAX_PATH (260 chars) — file extension truncated**

Pathway names are trimmed to 60 characters in saved filenames (plot titles
are unaffected). If you still see truncation, shorten `base_dir`.

---

## Differences from MATLAB v0.2

| Aspect | MATLAB v0.2 | R v0.2 |
|--------|-------------|--------|
| Input | `SingleCellExperiment` (scGEAtoolbox) | Seurat v4/v5 **or** SCE — auto-detected |
| NLS solver | Trust-Region (`fit()`) | Levenberg-Marquardt (`minpack.lm::nlsLM`) |
| Parallelism | `parfor` | `future.apply::future_lapply()` |
| Heatmap | `imagesc` | `pheatmap` |
| Gene plots | MATLAB `plot`/`violinplot` | `ggplot2` |
| GUI | `uicontrol` figures | Shiny web app (`app.R`) |
| Pathway pipeline | Not included | AUCell + pathway cosinor (Approaches A & B) |
| BH denominator | All genes (incl. zero-expr.) | Matched from v0.2 |
| Numerical agreement | — | < 0.1 % vs MATLAB (96.8 % gene overlap validated) |

---

## Version

**v0.2** — May 2026
