# TimeSCape v0.2 — Function Reference

Circadian rhythm detection in single-cell RNA-seq data.  
Three parallel implementations: **R**, **Python**, and **MATLAB**.  
Statistical results are numerically equivalent across all three platforms (< 0.1% difference).

---

## Repository Structure

```
circadian_ident_stat_test/
├── TimeSCape_R/                  # R implementation (Seurat / SingleCellExperiment)
│   ├── R/
│   │   ├── estimate_phaseR.R     # Core cosinor fit + significance test
│   │   ├── run_timescape.R       # Main pipeline engine
│   │   ├── generate_heatmap.R    # Heatmap + polar acrophase plot
│   │   ├── plot_gene.R           # Per-gene ZT expression plots
│   │   └── pathway_circadian.R   # Pathway enrichment, AUCell, GRN
│   ├── run_timescape_test.R      # Interactive single cell-type script
│   ├── run_full_pipeline_early.R
│   ├── run_full_pipeline_intermediate.R
│   └── run_full_pipeline_advanced.R
│
├── TimeSCape_py/                 # Python implementation (AnnData / Scanpy)
│   └── timescape/
│       ├── core.py               # Core cosinor fit + significance test
│       ├── pipeline.py           # Main pipeline engine
│       ├── normalize.py          # Library-size normalisation
│       ├── utils.py              # BH correction, ZT parsing, helpers
│       ├── visualize.py          # Heatmap + per-gene plots
│       ├── pathway.py            # Pathway enrichment (ORA + AUCell + cosinor)
│       └── grn.py                # Hub gene selection + GRN time series
│
└── TimeSCape_Matlab/             # MATLAB implementation (scGEAtoolbox SCE)
    ├── estimate_phaseR.m         # Core cosinor fit + significance test
    ├── sce_circ_phase_estimation_stattest.m  # Main pipeline engine
    ├── generateHeatmap_circ_simple.m
    ├── sce_circ_plot.m
    ├── sce_circ_plot_gene.m
    ├── bh_adjust_pvalues.m
    ├── progressbar.m
    └── TimeSCape_GUI.m           # MATLAB App Designer GUI
```

---

## R — Function Reference (`TimeSCape_R/R/`)

### `estimate_phaseR.R`

#### `estimate_phaseR(Xg_zts, actual_times, period12, test_type)`
Core cosinor fitting function. Fits a cosine model to single-cell expression data and tests significance via F-test or LRT.

| Argument | Type | Description |
|---|---|---|
| `Xg_zts` | named list of numeric vectors | Each element = expression values for all cells at one ZT time point. Empty elements (missing ZTs) are skipped automatically. |
| `actual_times` | numeric vector | True ZT hours for each slot in `Xg_zts`, e.g. `c(0, 3, 6, 9, 15, 18, 21)` when ZT12 is absent. |
| `period12` | logical | `TRUE` = 12-hr ultradian; `FALSE` = 24-hr circadian (default). |
| `test_type` | character | `"Ftest"` (default) or `"LRT"`. |

Returns a list: `acrophase`, `amp`, `mesor`, `period`, `p_value`, `rho`, `p_value_macro`, `R0` (per-ZT means). `amp` is always non-negative (canonical form); if the NLS converges with a negative amplitude, the sign is flipped and `acrophase` is shifted by `period/2` so it reflects the true peak time. `acrophase` is always wrapped into `[0, period)` via modulo.

---

### `run_timescape.R`

#### `build_tmeta(obj, zt_col)`
Parses ZT strings from Seurat or SingleCellExperiment metadata into numeric hours. Recognises formats: `"ZT00"`, `"ZT_03"`, `"ZT 6"`, plain `"0"`, `"3"`, etc.

| Argument | Type | Description |
|---|---|---|
| `obj` | Seurat or SCE | Input single-cell object. |
| `zt_col` | character | Metadata column containing ZT strings. |

Returns a `data.frame` with columns `zt_str` and `ZT_times`. Inspect and adjust `ZT_times` manually if automatic parsing fails.

```r
tmeta <- build_tmeta(seurat_obj, zt_col = "ZT_str")
```

---

#### `run_timescape(obj, celltype_col, zt_col, out_dir, period12, n_workers, ...)`
Main pipeline engine. Loops over all cell types in `celltype_col`, fits the cosinor model per gene, applies the dual significance criterion (F-test AND Pearson correlation), BH-corrects p-values, and writes output CSVs and optional heatmaps.

| Argument | Type | Default | Description |
|---|---|---|---|
| `obj` | Seurat or SCE | — | Input object. |
| `celltype_col` | character | — | Metadata column for cell type labels. |
| `zt_col` | character | — | Metadata column for ZT strings. |
| `out_dir` | character | — | Root output directory. |
| `period12` | logical | `FALSE` | 12-hr (`TRUE`) or 24-hr (`FALSE`) mode. |
| `n_workers` | integer | `1` | Parallel workers (`future.apply`). |
| `min_cells_per_zt` | integer | `5` | Minimum cells per ZT for a cell type to be analysed. |
| `plot_heatmap` | logical | `TRUE` | Generate heatmaps after fitting. |
| `use_normalized` | logical | `FALSE` | Use Seurat's pre-normalised `data` slot instead of raw counts. |

Returns a named list of result data frames (one per cell type). Also writes per-cell-type CSVs to `out_dir/{CellType}/01_gene_circadian/`.

```r
library(future)
plan(multisession, workers = 4)

results <- run_timescape(
  obj          = seurat_obj,
  celltype_col = "cell_type",
  zt_col       = "ZT_str",
  out_dir      = "TimeSCape_output_early",
  period12     = FALSE,
  n_workers    = 4
)
```

---

### `generate_heatmap.R`

#### `generate_heatmap(celltype, outdir, strict, period12, circ, custom_name, return_obj)`
Reads CSVs written by `run_timescape()` and renders a Z-score heatmap (genes × ZT time points), sorted by acrophase.

| Argument | Type | Default | Description |
|---|---|---|---|
| `celltype` | character | — | Sanitised cell type name (directory prefix). |
| `outdir` | character | — | Directory containing the cell-type subdirectory. |
| `strict` | logical | `TRUE` | `TRUE` = both F-test AND corr p < 0.05; `FALSE` = F-test only. |
| `period12` | logical | `FALSE` | Must match the run_timescape setting. |
| `circ` | logical | `FALSE` | Restrict to classical clock gene prefixes. |
| `custom_name` | character | `NULL` | Optional suffix for the PNG filename. |
| `return_obj` | logical | `FALSE` | Return the pheatmap object (for Shiny embedding). |

```r
generate_heatmap(
  celltype = "Tumor_cells",
  outdir   = "TimeSCape_output_early/Tumor_cells",
  strict   = TRUE
)
```

---

#### `load_stage_results(out_dir, period12, skip_failed)`
Reloads all CSVs from a completed pipeline run into a named list, without re-running the analysis. Use to regenerate or customise the polar acrophase plot.

```r
results <- load_stage_results("TimeSCape_output_early")
plot_clock_acrophase(results, stage = "Early")
```

---

#### `plot_clock_acrophase(results, stage, gene_list, cell_group_rules, strict)`
Polar (clock-face) chart of core clock gene acrophases across all cell types.

| Argument | Type | Default | Description |
|---|---|---|---|
| `results` | named list | — | From `run_timescape()` or `load_stage_results()`. |
| `stage` | character | — | Label for the plot title. |
| `gene_list` | character vector | `NULL` | `NULL` = built-in 11-gene clock set. |
| `cell_group_rules` | named list | `NULL` | `NULL` = built-in immune/tumour/structural bins. Custom: `list("Immune" = c("T cells", "B cells"))`. |
| `strict` | logical | `TRUE` | Only plot genes passing dual significance filter. |

---

### `plot_gene.R`

#### `plot_gene(tmeta, cust_cells, cust_gene, outdir, period12, print_scdata, sce, celltype_col, zt_col, use_violin)`
Returns a ggplot2 object showing the fitted cosine, per-ZT means, acrophase marker, and optionally single-cell overlay.

| Argument | Type | Default | Description |
|---|---|---|---|
| `tmeta` | data.frame | — | From `build_tmeta()`. |
| `cust_cells` | character | — | Sanitised cell type name. |
| `cust_gene` | character | — | Gene symbol. |
| `outdir` | character | — | Cell-type output directory. |
| `print_scdata` | logical | `FALSE` | Overlay individual cell expression. |
| `sce` | Seurat | `NULL` | Required when `print_scdata = TRUE`. |
| `use_violin` | logical | `FALSE` | Violin (`TRUE`) or dot scatter (`FALSE`). |

```r
p <- plot_gene(tmeta, "Tumor_cells", "Arntl",
               outdir = "TimeSCape_output_early/Tumor_cells")
print(p)
```

---

### `pathway_circadian.R`

#### `choose_collection(default_collection, default_subcategory)`
Interactively browse MSigDB collections and pick one. In non-interactive sessions, returns the supplied defaults.

```r
coll <- choose_collection()  # interactive
# or
coll <- list(collection = "C2", subcategory = "CP:KEGG_LEGACY")
```

---

#### `pull_genesets(collection, subcategory, species, min_size, max_size, deduplicate)`
Downloads gene sets from MSigDB via `msigdbr` and removes redundant parent/child terms.

| Argument | Type | Default | Description |
|---|---|---|---|
| `collection` | character | `"C2"` | MSigDB collection: `"C2"` (curated), `"C5"` (GO). |
| `subcategory` | character | — | E.g. `"CP:KEGG_LEGACY"`, `"CP:REACTOME"`, `"GO:BP"`. |
| `species` | character | `"Mus musculus"` | Species name as in msigdbr. |
| `min_size` | integer | `10` | Minimum gene set size. |
| `max_size` | integer | `500` | Maximum gene set size. |
| `deduplicate` | logical | `TRUE` | Remove positive/negative regulation variants. |

```r
genesets <- pull_genesets("C2", "CP:KEGG_LEGACY")
```

---

#### `inspect_genesets(gs_list, outfile)`
Diagnostic: shows all gene set names, sizes, and deduplication status. Useful before running ORA.

```r
inspect_genesets(genesets, outfile = "genesets_summary.csv")
```

---

#### `phase_bin_analysis(conf_genes, genesets, bin_width, n_top, min_overlap, min_bin_genes, p_thresh, use_padj, exclude_patterns)`
Bins confident circadian genes by acrophase into ZT windows, runs ORA per bin, and builds phase-restricted gene sets for AUCell.

| Argument | Type | Default | Description |
|---|---|---|---|
| `conf_genes` | data.frame | — | Confident gene table from `run_timescape()`. |
| `genesets` | named list | — | From `pull_genesets()`. |
| `bin_width` | numeric | `2` | ZT window width in hours. |
| `n_top` | integer | `5` | Top pathways returned per bin. |
| `min_overlap` | integer | `3` | Minimum gene overlap for a pathway hit. |
| `min_bin_genes` | integer | `5` | Minimum genes in a bin to attempt ORA. |
| `p_thresh` | numeric | `0.05` | Significance threshold. |
| `use_padj` | logical | `TRUE` | Use BH-adjusted p-value. |
| `exclude_patterns` | character vector | `c("DISEASE","VIRAL","INFECTION")` | Filter out pathway names matching these strings. |

Returns `list(bin_table, ora_results, phase_gs)`.

```r
pb <- phase_bin_analysis(conf_genes = results[["Tumor_cells"]]$confident,
                          genesets   = genesets)
```

---

#### `auc_score_cells(sce, gene_sets, celltype_col, celltype, zt_col, ...)`
Runs AUCell on a named gene-set list for a specific cell type, returning a matrix of per-cell AUC scores.

---

#### `pathway_cosinor(auc_matrix, tmeta, period12, fdr_thresh)`
Applies the TimeSCape cosinor framework to pathway AUCell scores. Each pathway becomes the "gene" and each cell's AUC score is its expression. Returns a table of significantly oscillating pathways with acrophase, amplitude, and p-values.

```r
pathway_results <- pathway_cosinor(auc_matrix, tmeta, period12 = FALSE)
```

---

#### `select_hub_genes(sce, conf_genes, clock_genes_ref, celltype_col, celltype, cor_thresh, p_thresh, hub_pct, min_hub)`
Identifies hub genes for GRN construction using global Pearson correlation across all ZT time points (pooled).

| Argument | Type | Default | Description |
|---|---|---|---|
| `conf_genes` | character vector | — | Confident circadian gene symbols. |
| `clock_genes_ref` | character vector | — | Core clock genes to always include. |
| `cor_thresh` | numeric | `0.30` | Minimum \|r\| for an edge to count toward degree. |
| `hub_pct` | numeric | `0.10` | Top fraction of genes by degree (0.10 = top 10%). |
| `min_hub` | integer | `5` | Minimum hubs returned; threshold auto-relaxed if needed. |

Returns `list(genes, cor_mat, adj_mat, degree, hub_genes, n_cells)`.

---

#### `plot_grn_timeseries(sce, hub_genes, pathway_genes, tmeta, celltype_col, celltype, zt_col, title, outfile)`
Builds one co-expression network panel per ZT time point using only cells from that slice. Node layout is fixed from the pooled (all-ZT) network for visual stability across panels.

---

#### `write_pathway_results(results, outpath, celltype)`
Writes pathway cosinor results to Excel.

---

## Python — Function Reference (`TimeSCape_py/timescape/`)

### `core.py`

#### `estimate_phase_r(Xg_zts, actual_times, period, test_type)`
Core cosinor fitting function. Python equivalent of `estimate_phaseR.R` and `estimate_phaseR.m`.

| Argument | Type | Default | Description |
|---|---|---|---|
| `Xg_zts` | list of np.ndarray | — | Per-ZT expression arrays. Empty arrays skipped automatically. |
| `actual_times` | np.ndarray | — | Numeric ZT hours for each slot. |
| `period` | float | `24.0` | Period in hours (24 or 12). |
| `test_type` | str | `"Ftest"` | `"Ftest"` or `"LRT"`. |

Returns dict: `acrophase`, `amp`, `mesor`, `period`, `p_value`, `rho`, `p_value_macro`, `R0`. `amp` is always non-negative (canonical form); `acrophase` is always in `[0, period)` via modulo wrapping.

```python
from timescape.core import estimate_phase_r
result = estimate_phase_r(Xg_zts, actual_times=np.array([0,3,6,9,12,15,18,21]))
```

---

### `pipeline.py`

#### `run_timescape(adata, celltype_col, zt_col, period, test_type, fdr_thresh, min_cells_per_zt)`
Main pipeline engine. Loops over cell types in `adata.obs[celltype_col]`, fits cosinor per gene, applies dual significance filter, and BH-corrects p-values.

| Argument | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | — | Scanpy AnnData object with raw counts in `.X`. |
| `celltype_col` | str | — | `obs` column for cell type labels. |
| `zt_col` | str | — | `obs` column for ZT strings. |
| `period` | float | `24.0` | Circadian period. |
| `test_type` | str | `"Ftest"` | `"Ftest"` or `"LRT"`. |
| `fdr_thresh` | float | `0.05` | BH-corrected p-value threshold for confident genes. |
| `min_cells_per_zt` | int | `5` | Minimum cells per ZT to include a time point. |

Returns a dict of DataFrames keyed by cell type.

```python
from timescape.pipeline import run_timescape
results = run_timescape(adata, celltype_col="cell_type", zt_col="ZT_str")
```

---

### `normalize.py`

#### `normalize_lib_size(X, target)`
Normalises a genes × cells count matrix to `target` counts per cell and applies log1p.

#### `normalize_adata(adata, target, layer)`
In-place library-size normalisation of an AnnData object. Stores result in `adata.layers[layer]`.

---

### `utils.py`

#### `bh_adjust(pvalues)`
Benjamini-Hochberg FDR correction. Returns adjusted p-values as `np.ndarray`.

#### `wrap_acrophase(acro, period)`
Wraps an acrophase value to `[0, period)`.

#### `parse_zt_string(s)`
Parses a ZT string (`"ZT03"`, `"ZT_6"`, `"zt12"`, `"3"`) to a numeric float. Returns `None` on failure.

#### `build_tmeta(zt_strings, exclude)`
Builds a ZT metadata DataFrame from a list of ZT strings. Returns a `pd.DataFrame` with columns `zt_str` and `ZT_times`, sorted by `ZT_times`.

```python
from timescape.utils import build_tmeta
tmeta = build_tmeta(adata.obs["ZT_str"].unique().tolist())
```

---

### `visualize.py`

#### `generate_heatmap(result_df, tmeta, celltype, outdir, strict, period12)`
Blue-white-red Z-score heatmap of circadian genes × ZT time points, sorted by acrophase.

#### `plot_gene_single(result_df, adata, gene, celltype, tmeta, celltype_col, zt_col, outdir, show_cells, use_violin)`
Per-gene plot: fitted cosine + per-ZT means + optional single-cell overlay.

#### `save_batch_plots(result_df, adata, celltype, tmeta, celltype_col, zt_col, outdir, n_top)`
Saves per-gene plots for the top `n_top` confident genes.

---

### `pathway.py`

Python equivalent of `pathway_circadian.R`. Requires: `pip install gseapy decoupler openpyxl`.

#### `pull_genesets(organism, collection, min_size, max_size)`
Fetches gene sets from Enrichr/MSigDB via gseapy. Mirrors `pull_genesets()` in R.

| Argument | Default | Description |
|---|---|---|
| `organism` | `"mouse"` | `"mouse"` or `"human"`. |
| `collection` | `"KEGG_2019_Mouse"` | Any Enrichr library name (see `gseapy.get_library_name()`). Examples: `"KEGG_2019_Mouse"`, `"Reactome_2022"`, `"GO_Biological_Process_2021"` for mouse; `"KEGG_2021_Human"` for human. |
| `min_size` | `10` | Minimum gene set size after filtering. |
| `max_size` | `500` | Maximum gene set size. |

Returns `dict` mapping pathway name → list of gene symbols.

```python
from timescape.pathway import pull_genesets
genesets = pull_genesets("mouse", "KEGG_2019_Mouse")
```

---

#### `phase_bin_analysis(conf_df, genesets, universe, ...)`
Bins confident circadian genes by acrophase into ZT windows and runs ORA (hypergeometric test, equivalent to `clusterProfiler::enricher()`) per bin. Returns phase-restricted gene sets for AUCell (Method A).

| Argument | Default | Description |
|---|---|---|
| `conf_df` | — | Confident gene DataFrame from `run_timescape()`. |
| `genesets` | — | From `pull_genesets()`. |
| `universe` | — | List of all genes tested (not genome-wide — avoids p-value inflation). |
| `bin_width` | `2.0` | ZT window in hours. |
| `n_top` | `5` | Top pathway hits retained per bin. |
| `min_overlap` | `3` | Minimum gene overlap. |
| `p_thresh` | `0.05` | Significance threshold. |
| `use_padj` | `True` | Use BH-adjusted p-value. |
| `exclude_patterns` | `["DISEASE","VIRAL","INFECTION"]` | Pathway name filters. |

Returns `dict(bin_table, ora_results, phase_gs)`.

```python
from timescape.pathway import phase_bin_analysis
pb = phase_bin_analysis(conf_df, genesets, universe=tested_genes)
# pb["phase_gs"] → {"ZT04.0-06.0__KEGG_CIRCADIAN_RHYTHM": ["Arntl", "Clock", ...], ...}
```

---

#### `auc_score_cells(adata, genesets, min_gs_size, auc_max_rank)`
AUCell pathway activity scoring (pathways × cells). Uses `decoupler.run_aucell()` if available; falls back to a numpy implementation.

| Argument | Default | Description |
|---|---|---|
| `adata` | — | Log-normalised AnnData (cells × genes). |
| `genesets` | — | From `pull_genesets()` or `phase_bin_analysis()["phase_gs"]`. |
| `min_gs_size` | `5` | Minimum pathway genes overlapping the matrix. |
| `auc_max_rank` | `0.05` | Top fraction of genes used as AUC threshold (5%). |

Returns `np.ndarray` (pathways × cells) with `.pathway_names` attribute.

```python
from timescape.pathway import auc_score_cells
auc_A = auc_score_cells(adata, pb["phase_gs"])     # Method A (phase-restricted)
auc_B = auc_score_cells(adata, genesets)            # Method B (full pathways)
```

---

#### `pathway_cosinor(auc_mat, pathway_names, obs, tmeta, target_ct, ...)`
Applies the same cosinor framework as `run_timescape()` to pathway AUCell scores. Each pathway's per-cell AUC score is treated as its expression.

Returns `dict(stats, zt_means)` where `stats` has the same columns as the T1 gene table.

```python
from timescape.pathway import pathway_cosinor
res = pathway_cosinor(auc_A, auc_A.pathway_names, adata.obs, tmeta,
                      target_ct="Tumor_cells")
```

---

#### `write_pathway_results(results, outpath, celltype)`
Saves pathway cosinor stats to a formatted Excel workbook (two sheets: All_Pathways and Confident_Pathways, with header colour and confident-row highlighting).

```python
from timescape.pathway import write_pathway_results
write_pathway_results(res, "Tumor_cells_pathway_A.xlsx", celltype="Tumor_cells")
```

---

### `grn.py`

Python equivalent of the GRN section in `pathway_circadian.R`. Requires: `pip install networkx matplotlib`.

#### `select_hub_genes(X_ct, gene_names, cor_thresh, p_thresh, hub_pct, min_hub)`
Selects hub genes from a circadian gene pool using global Pearson correlation on all-ZT pooled cells.

| Argument | Default | Description |
|---|---|---|
| `X_ct` | — | Dense expression matrix (cells × genes or genes × cells) for one cell type. |
| `gene_names` | — | Gene labels (columns or rows of `X_ct`). |
| `cor_thresh` | `0.30` | Minimum \|r\| for an edge to count toward degree. |
| `p_thresh` | `0.05` | Maximum p-value for an edge. |
| `hub_pct` | `0.10` | Top fraction of genes by degree (0.10 = top 10%). |
| `min_hub` | `5` | Minimum hubs; threshold is relaxed automatically if needed. |

Returns `dict(genes, cor_mat, adj_mat, degree, hub_genes, n_cells)`.

```python
from timescape.grn import select_hub_genes
hub_info = select_hub_genes(X_ct, gene_names=conf_genes + clock_genes)
hub_genes = hub_info["hub_genes"]
```

---

#### `plot_grn_timeseries(X_gc, gene_names, zt_v, hub_genes, actual_times, ...)`
Builds one co-expression network panel per ZT time point (networkx + matplotlib). Node layout is fixed from the pooled network for visual stability across panels. Edge colour encodes signed Pearson r (red = positive, blue = negative).

| Argument | Default | Description |
|---|---|---|
| `X_gc` | — | Genes × cells matrix for the target cell type, all ZT pooled. |
| `gene_names` | — | Gene labels for rows of `X_gc`. |
| `zt_v` | — | Per-cell ZT hours. |
| `hub_genes` | — | Genes to display (3–15 for readable plots). |
| `actual_times` | — | Unique ZT hours to plot (one panel each). |
| `cor_thresh` | `0.20` | Edge threshold for drawn edges. |
| `outpath` | `None` | Save figure to this path (PNG/PDF) if given. |

Returns a `matplotlib.figure.Figure`.

```python
from timescape.grn import plot_grn_timeseries
fig = plot_grn_timeseries(X_gc, gene_names, zt_v, hub_genes,
                           actual_times=np.array([0,3,6,9,12,15,18,21]),
                           title="Tumor cells — GRN", outpath="grn_tumor.png")
```

---

### Full Python pathway + GRN example

```python
import scanpy as sc
import numpy as np
from timescape import (run_timescape, build_tmeta,
                        pull_genesets, phase_bin_analysis,
                        auc_score_cells, pathway_cosinor, write_pathway_results,
                        select_hub_genes, plot_grn_timeseries)

adata  = sc.read_h5ad("my_data.h5ad")
tmeta  = build_tmeta(adata.obs["ZT_str"].unique().tolist())

# Step A: gene-level circadian
T1, T2 = run_timescape(adata, tmeta, celltype_col="cell_type", zt_col="ZT_str")
conf    = T1[T1.cell_type == "Tumor_cells"]
conf_df = conf[(conf.pvalue < 0.05) & (conf.pvalue_corr < 0.05)]

# Step B: pathway enrichment (Method A — phase-restricted ORA + AUCell)
genesets  = pull_genesets("mouse", "KEGG_2019_Mouse")
pb        = phase_bin_analysis(conf_df, genesets,
                                universe=T1[T1.cell_type=="Tumor_cells"]["Genes"].tolist())
auc_A     = auc_score_cells(adata, pb["phase_gs"])
res_A     = pathway_cosinor(auc_A, auc_A.pathway_names, adata.obs, tmeta,
                             target_ct="Tumor_cells")
write_pathway_results(res_A, "Tumor_cells_approachA.xlsx")

# Step C: Method B — full pathway AUCell
auc_B  = auc_score_cells(adata, genesets)
res_B  = pathway_cosinor(auc_B, auc_B.pathway_names, adata.obs, tmeta,
                          target_ct="Tumor_cells")
write_pathway_results(res_B, "Tumor_cells_approachB.xlsx")

# Step D: GRN
ct_mask   = adata.obs["cell_type"] == "Tumor_cells"
X_ct      = adata[ct_mask].X.toarray().T          # genes × cells
zt_v      = np.array([tmeta.set_index("old_labels").loc[z, "ZT_times"]
                       for z in adata.obs.loc[ct_mask, "ZT_str"]])
clock_ref = ["Arntl", "Clock", "Per1", "Per2", "Cry1", "Cry2",
             "Nr1d1", "Nr1d2", "Dbp", "Rora", "Bhlhe40"]
gene_pool = conf_df["Genes"].tolist() + [g for g in clock_ref if g in adata.var_names]
hub_info  = select_hub_genes(X_ct.T, gene_names=list(adata.var_names))  # pass full, filter internally
fig = plot_grn_timeseries(X_ct, list(adata.var_names), zt_v,
                           hub_info["hub_genes"],
                           actual_times=np.sort(np.unique(zt_v)),
                           title="Tumor cells GRN", outpath="grn_tumor.png")
```

---

## MATLAB — Function Reference (`TimeSCape_Matlab/`)

### `estimate_phaseR.m`

#### `[acrophase, amp, period, mesor, p_value, rho, p_value_macro] = estimate_phaseR(Xg_zts, actual_times, period12, test_type)`
Core cosinor fitting function. Accepts missing time points via `actual_times`.

| Argument | Description |
|---|---|
| `Xg_zts` | Cell array (1 × nzts). Each cell = expression values at one ZT. Empty cells skipped. |
| `actual_times` | Numeric vector of true ZT hours, e.g. `[0 3 6 9 15 18 21]`. |
| `period12` | Logical; `true` = 12-hr period. |
| `test_type` | `'Ftest'` or `'LRT'`. |

`amp` is always non-negative in the output (canonical form); `acrophase` is always in `[0, period)` via modulo.

---

### `sce_circ_phase_estimation_stattest.m`

#### `[T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, rm_low_conf, period12, custom_genelist, custom_celltype, plot_heat, norm_str, num_cores)`
Main pipeline engine. Loops over all cell types, fits cosinor per gene, applies dual significance filter, BH-corrects, and writes CSVs.

| Argument | Default | Description |
|---|---|---|
| `sce` | — | SingleCellExperiment object. |
| `tmeta` | — | Table: `old_labels \| new_labels \| ZT_times`. |
| `rm_low_conf` | `true` | Write confident-only secondary CSVs. |
| `period12` | `false` | 12-hr ultradian mode. |
| `custom_genelist` | `[]` | Restrict to specific genes. |
| `custom_celltype` | `[]` | Restrict to one cell type. |
| `plot_heat` | `true` | Generate heatmap after analysis. |
| `norm_str` | `'lib_size'` | `'lib_size'` \| `'none'` \| `'magic_impute'`. |
| `num_cores` | auto | Parallel pool workers. |

```matlab
times      = [0 3 6 9 12 15 18 21]';
old_labels = unique(sce.c_batch_id);
new_labels = string(arrayfun(@(t) sprintf('ZT%02d',t), times, 'UniformOutput', false));
tmeta      = table(old_labels, new_labels, times, 'VariableNames', {'old_labels','new_labels','ZT_times'});
[T1, T2]   = sce_circ_phase_estimation_stattest(sce, tmeta);
```

Outputs (per cell type):
- `*_circadian_analysis_all.csv` — all tested genes
- `*_circadian_analysis_confident.csv` — dual-test significant genes
- `*_circadian_ZTs_mean.csv` — raw per-ZT means
- `*_circadian_ZTs_mean_normalized.csv` — ZT00-normalised means

---

### `generateHeatmap_circ_simple.m`

#### `generateHeatmap_circ_simple(celltype, outdir, strict, period12)`
Blue-white-red Z-score heatmap of confident circadian genes, sorted by acrophase.

---

### `sce_circ_plot.m` / `sce_circ_plot_gene.m`

Per-cell-type summary plots and individual gene expression plots with fitted cosine overlay.

---

### `bh_adjust_pvalues.m`

#### `padj = bh_adjust_pvalues(pvals)`
Benjamini-Hochberg FDR correction. Equivalent to R's `p.adjust(method="BH")`.

---

### `TimeSCape_GUI.m`

MATLAB App Designer GUI for interactive analysis. Launch with:
```matlab
TimeSCape_GUI
```
Supports light/dark mode. Embeds heatmap and gene plots inline.

---

## Quick-Start Examples

### R — Run full pipeline
```r
library(Seurat)
library(future)
library(future.apply)
source("TimeSCape_R/R/estimate_phaseR.R")
source("TimeSCape_R/R/run_timescape.R")
source("TimeSCape_R/R/generate_heatmap.R")

plan(multisession, workers = 4)

seurat_obj <- readRDS("my_data.rds")

results <- run_timescape(
  obj          = seurat_obj,
  celltype_col = "cell_type",
  zt_col       = "ZT_str",
  out_dir      = "TimeSCape_output_early",
  period12     = FALSE,
  n_workers    = 4
)
```

### Python — Run pipeline
```python
import scanpy as sc
from timescape.pipeline import run_timescape

adata = sc.read_h5ad("my_data.h5ad")
results = run_timescape(adata, celltype_col="cell_type", zt_col="ZT_str")
```

### MATLAB — Run pipeline
```matlab
load('my_sce.mat');  % loads sce
times      = [0 3 6 9 12 15 18 21]';
old_labels = unique(sce.c_batch_id);
new_labels = string(arrayfun(@(t) sprintf('ZT%02d',t), times, 'UniformOutput', false));
tmeta      = table(old_labels, new_labels, times, 'VariableNames', {'old_labels','new_labels','ZT_times'});
[T1, T2]   = sce_circ_phase_estimation_stattest(sce, tmeta);
```

---

## Dependencies

### R
```r
install.packages(c("Seurat", "future", "future.apply", "minpack.lm",
                   "ggplot2", "gridExtra", "openxlsx", "dplyr",
                   "igraph", "ggraph", "cowplot", "msigdbr"))
BiocManager::install(c("clusterProfiler", "AUCell", "SingleCellExperiment"))
```

### Python
```bash
conda env create -f TimeSCape_py/environment.yml
conda activate timescape

# Additional dependencies for pathway enrichment and GRN:
pip install gseapy decoupler openpyxl networkx matplotlib
```

Core statistical pipeline (`run_timescape`, `estimate_phase_r`) requires only numpy, scipy, pandas, joblib, and scanpy. The pathway and GRN modules (`pathway.py`, `grn.py`) require the packages above and are imported lazily — the core pipeline works without them.

### MATLAB
- MATLAB R2023a or later
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox
- Curve Fitting Toolbox
- [scGEAtoolbox](https://github.com/jamesjcai/scGEAToolbox)

---

## Contact

Selim Romero — ssromerogon@tamu.edu
