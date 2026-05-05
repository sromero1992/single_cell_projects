# TimeSCape R — Function Reference

Source files in this folder. See `../README.md` for installation,
setup, and full usage examples.

---

## Files

| File | Purpose |
|------|---------|
| `estimate_phaseR.R` | Core cosinor fitting for a single gene |
| `run_timescape.R` | Gene-level pipeline — loops over cell types |
| `generate_heatmap.R` | Reads CSVs and renders the z-score heatmap |
| `plot_gene.R` | Single-gene plot (returns ggplot object) |
| `pathway_circadian.R` | AUCell scoring + pathway cosinor + plots |

---

## `estimate_phaseR.R`

### `estimate_phaseR(Xg_zts, actual_times, period12, test_type)`

Fits `A·cos(2π(t − φ)/T) + M` to single-cell expression data at multiple
ZT time points and tests significance.

| Argument | Type | Description |
|----------|------|-------------|
| `Xg_zts` | list of numeric vectors | One element per ZT; empty elements skipped |
| `actual_times` | numeric vector | True ZT hours matching each slot |
| `period12` | logical | `TRUE` = 12-hr; `FALSE` = 24-hr |
| `test_type` | `"Ftest"` / `"LRT"` | Use `"Ftest"` (matches MATLAB) |

Returns a named list: `acrophase`, `amp`, `period`, `mesor`, `p_value`,
`rho`, `p_value_macro`. Returns all `NA` on fit failure.

### `diagnose_phaseR(obj, celltype_col, zt_col, tmeta, target_ct, norm_str, n)`

Runs `estimate_phaseR` on `n` random genes and prints the raw output.
Use to debug "0 genes tested" or all-NA results.

---

## `run_timescape.R`

### `build_tmeta(obj, zt_col)` / `build_tmeta_from_seurat(obj, zt_col)`

Parses unique ZT strings from a metadata column into numeric hours.
Returns a `data.frame` with columns `zt_str` and `ZT_times`, sorted by hour.
`build_tmeta_from_seurat` is a backward-compatible alias.

### `run_timescape(obj, celltype_col, zt_col, tmeta, ...)`

Main gene-level pipeline. For each (cell type × group) combination:
- Extracts and optionally re-normalises the expression submatrix
- Fits cosinor model to every gene in parallel
- Applies BH correction separately to F-test and correlation p-values
- Writes six CSVs + optional heatmap to `outdir/{combo_name}/`
- Returns a named list: `$T1` (stats) and `$T2` (per-ZT means)

**Key arguments**

| Argument | Default | Notes |
|----------|---------|-------|
| `celltype_col` | `"cell_type"` | Metadata column for cell-type labels |
| `zt_col` | `"ZT_str"` | Metadata column for ZT strings |
| `tmeta` | `NULL` | Built automatically if NULL |
| `norm_str` | `"lib_size"` | `"lib_size"` · `"logcounts"` · `"none"` |
| `period12` | `FALSE` | 12-hr or 24-hr period |
| `rm_low_conf` | `TRUE` | Write confident-only CSV subsets |
| `plot_heat` | `FALSE` | Save heatmap PNG after each cell type |
| `custom_celltype` | `NULL` | Restrict to specific cell types |
| `group_col` | `NULL` | Optional 2nd split variable |
| `n_workers` | `1L` | Parallel workers |
| `outdir` | `getwd()` | Root output folder |

---

## `generate_heatmap.R`

### `generate_heatmap(celltype, strict, custom_name, circ, period12, outdir, return_obj)`

Reads `*_circadian_analysis_all.csv` and `*_circadian_ZTs_mean.csv`,
filters genes, row-normalises (z-score), and renders a blue→white→red
heatmap via `pheatmap`. Genes sorted by acrophase then amplitude.
Negative-amplitude rows are flipped before sorting.

| Argument | Default | Notes |
|----------|---------|-------|
| `celltype` | — | Sanitised combo name matching directory/file prefix |
| `strict` | `TRUE` | Keep only genes with F-test AND corr p < 0.05 |
| `custom_name` | `""` | Optional suffix appended to output filenames |
| `circ` | `FALSE` | Restrict to classical circadian gene prefixes |
| `period12` | `FALSE` | Read 12-hr or 24-hr files |
| `outdir` | `getwd()` | Cell-type combo subdirectory (where CSVs live) |
| `return_obj` | `FALSE` | If `TRUE`, return pheatmap object for Shiny rendering |

---

## `plot_gene.R`

### `plot_gene_single(tmeta, cust_cells, period12, cust_gene, ...)`

Reads pre-computed CSVs and returns a `ggplot2` object showing:
- Blue cosine fit line
- Orange per-ZT mean markers
- Dashed acrophase line (labelled on x-axis)
- Optional single-cell overlay (violin or jitter)

Does **not** re-run the analysis. Requires `run_timescape()` to have been
run first (reads CSVs from `outdir`).

| Argument | Default | Notes |
|----------|---------|-------|
| `cust_cells` | — | Combo name matching directory/file prefix |
| `cust_gene` | — | Gene symbol string |
| `print_scdata` | `FALSE` | Overlay individual cell expression |
| `sce` | `NULL` | Seurat object — required when `print_scdata = TRUE` |
| `use_violin` | `FALSE` | `TRUE` = violin; `FALSE` = jitter dots |
| `outdir` | `getwd()` | Cell-type combo subdirectory |

Returns a `ggplot2` object — add layers or pass directly to `ggsave()`.

---

## `pathway_circadian.R`

### `pull_genesets(collection, subcategory, species, min_size, max_size, deduplicate)`

Downloads gene sets from MSigDB via `msigdbr` and returns a named list
(pathway → character vector of gene symbols), filtered by size.

### `auc_score_cells(obj, genesets, use_norm, auc_max_rank, n_cores, min_gs_size)`

Scores all cells for all gene sets using `AUCell`. Ranks genes per cell
once, then evaluates each gene set as a fast lookup. Returns a
**pathways × cells** numeric matrix. Results are deterministic and
cell-type-independent (per-cell ranking).

| Argument | Default | Notes |
|----------|---------|-------|
| `obj` | — | Seurat object |
| `genesets` | — | Named list from `pull_genesets()` or `msigdbr` |
| `use_norm` | `TRUE` | Use `data` slot (normalised); `FALSE` = `counts` |
| `auc_max_rank` | `0.05` | Top fraction of ranked genes used for AUC |
| `n_cores` | `1L` | AUCell parallel cores |
| `min_gs_size` | `5L` | Skip gene sets with fewer measured genes |

### `pathway_cosinor(auc_mat, meta, celltype_col, zt_col, tmeta, target_ct, period12)`

Fits the cosinor model to each row of `auc_mat` (pathway scores across
cells of `target_ct`), averaging scores per ZT time point.
Returns a list with `$stats` (one row per pathway, same columns as
gene-level stats table) and `$score_mat` (ZT-averaged scores).

> **Important**: `auc_mat` must have been scored on the **subsetted**
> cells only (same cell type as `target_ct`). `meta` must match those
> same cells. Passing the full-object metadata causes a subscript error.

### `write_pathway_results(path_results, xlsx_path, celltype)`

Writes two sheets to an Excel file: `All_Pathways` and
`Confident_Pathways` (p < 0.05 and pvalue_corr < 0.05).

### `plot_pathway_single(auc_mat, path_results, meta, celltype_col, zt_col, tmeta, target_ct, target_pathway, period12, use_violin)`

Returns a `ggplot2` cosine-fit plot for one pathway (same style as
`plot_gene_single`). Pass to `ggsave()` or include in a grid.

### `save_batch_pathway_plots(auc_mat, path_results, meta, ..., n_top, outdir)`

Saves one PNG per confident pathway (sorted by adjusted p-value).
Filenames are trimmed to 60 characters to avoid Windows MAX_PATH issues.
