# TimeSCape R — Function Reference

This folder contains the four source files that make up the R pipeline.
See `../README.md` for installation, setup, and full usage examples.

---

## Files

| File | Purpose |
|------|---------|
| `estimate_phaseR.R` | Core cosinor fitting for a single gene at one time |
| `run_timescape.R` | Main pipeline — loops over cell types (and groups) |
| `generate_heatmap.R` | Reads CSVs and renders the z-score heatmap |
| `plot_gene.R` | Single-gene plot + batch PNG export |

---

## `estimate_phaseR.R`

### `estimate_phaseR(Xg_zts, actual_times, period12, test_type)`

Fits `A·cos(2π(t − φ)/T) + M` to single-cell expression data at multiple ZT
time points and tests significance.

**Arguments**

| Argument | Type | Description |
|----------|------|-------------|
| `Xg_zts` | list of numeric vectors | One element per ZT time point, each holding all cell expression values at that time. Empty elements are skipped automatically. |
| `actual_times` | numeric vector | True ZT hours matching each slot of `Xg_zts`. Works with missing time points — no imputation. |
| `period12` | logical | `TRUE` = 12-hr period; `FALSE` = 24-hr. |
| `test_type` | `"Ftest"` or `"LRT"` | Significance test. Use `"Ftest"` (matches MATLAB). |

**Returns** a named list:

| Name | Description |
|------|-------------|
| `acrophase` | Estimated peak time (hrs, within [0, period]) |
| `amp` | Cosine amplitude (signed) |
| `period` | Period used (12 or 24) |
| `mesor` | Midline (mean expression level) |
| `p_value` | F-test p-value |
| `rho` | Pearson r between per-ZT means and fitted curve |
| `p_value_macro` | p-value for that Pearson correlation |

Returns all `NA` on fit failure (< 4 time points, or NLS non-convergence).

---

## `run_timescape.R`

### `build_tmeta_from_seurat(seurat_obj, zt_col)`

Parses unique ZT strings from a Seurat metadata column into numeric hours.
Supports `"ZT00"`, `"ZT03"`, `"zt12"`, `"ZT_06"`, `"ZT 3"`, plain `"0"` …

Returns a `data.frame` with columns `zt_str` and `ZT_times`, sorted by hour.
Inspect and edit `ZT_times` manually if auto-parsing fails.

---

### `run_timescape(seurat_obj, ...)`

Main analysis pipeline. Loops over cell types (and optionally a second group
column), fits the cosinor model to every gene, applies BH correction, and
writes 6 CSVs + optional heatmap per combination.

**Key arguments**

| Argument | Default | Notes |
|----------|---------|-------|
| `celltype_col` | `"cell_type"` | Metadata column for cell-type labels |
| `zt_col` | `"ZT_str"` | Metadata column for ZT strings |
| `tmeta` | `NULL` | Built automatically if NULL |
| `norm_str` | `"lib_size"` | `"lib_size"` · `"seurat"` · `"none"` |
| `period12` | `FALSE` | 12-hr or 24-hr period |
| `group_col` | `NULL` | Optional 2nd split variable (cancer stage, replicate …) |
| `custom_group` | `NULL` | Restrict to specific group values |
| `outdir` | `getwd()` | Root output folder |

**Normalization note** — `"lib_size"` normalises per-cell inside the loop so
peak RAM equals one dense (cell_type × group) submatrix, not the full object.
This is bit-for-bit identical to normalising the full matrix first because
library-size is a per-cell operation.

**Output directory structure**

```
outdir/
  CellType/                          # no group_col
  CellType_GroupValue/               # with group_col
    CellType_period_24_circadian_analysis_all.csv
    CellType_period_24_circadian_analysis_confident.csv
    CellType_period_24_circadian_ZTs_mean.csv
    CellType_period_24_circadian_ZTs_mean_normalized.csv
    CellType_period_24_circadian_ZTs_mean_confident.csv
    CellType_period_24_circadian_ZTs_mean_normalized_confident.csv
    CellType_period_24_heatmap_strict.png
  all_cell_types_period_24_summary_results.csv
```

**Returns** a named list keyed by `combo_name` (e.g. `"Hepatocytes"` or
`"Hepatocytes_early"`), each element holding `$T1` (stats table) and `$T2`
(per-ZT means).

---

## `generate_heatmap.R`

### `generate_heatmap(celltype, strict, custom_name, circ, period12, outdir, return_obj)`

Reads the `*_circadian_analysis_all.csv` and `*_circadian_ZTs_mean.csv` files
written by `run_timescape()`, filters genes, row-normalises (z-score), and
renders a blue→white→red heatmap via `pheatmap`.

| Argument | Default | Notes |
|----------|---------|-------|
| `celltype` | — | Sanitised combo name matching directory/file prefix |
| `strict` | `TRUE` | Keep only genes with F-test AND corr p < 0.05 |
| `custom_name` | `""` | Optional suffix appended to output filenames |
| `circ` | `FALSE` | Restrict to classical circadian gene prefixes |
| `period12` | `FALSE` | Read 12-hr or 24-hr files |
| `outdir` | `getwd()` | Cell-type combo subdirectory (where CSVs live) |
| `return_obj` | `FALSE` | If `TRUE`, return the pheatmap object for Shiny rendering |

Genes are sorted by acrophase then amplitude. Negative-amplitude rows are
flipped (amplitude → positive, acrophase shifted +12 h) before sorting.

---

## `plot_gene.R`

### `plot_gene_single(tmeta, cust_cells, period12, cust_gene, ...)`

Reads pre-computed CSVs and returns a `ggplot2` object showing:
- Blue cosine fit line
- Orange per-ZT mean ± markers
- Dark-red dashed acrophase line (labelled on axis)
- Optional single-cell overlay — violin or jitter scatter

Does **not** re-run the analysis. Run `run_timescape()` first.

| Argument | Default | Notes |
|----------|---------|-------|
| `cust_cells` | — | Combo name matching directory/file prefix |
| `cust_gene` | — | Gene name string |
| `print_scdata` | `FALSE` | Overlay individual cell expression |
| `sce` | `NULL` | Seurat object — required when `print_scdata = TRUE` |
| `use_violin` | `FALSE` | `TRUE` = violin; `FALSE` = jitter dots |
| `outdir` | `getwd()` | Cell-type combo subdirectory |

Returns a `ggplot2` object — add layers or pass to `ggsave()` freely.

---

### `save_batch_plots(tmeta, cust_cells, plot_type, period12, outdir)`

Loops over a filtered gene list and saves one PNG per gene.

| `plot_type` | Gene set |
|-------------|----------|
| `1` | Confident (F-test AND corr p < 0.05) |
| `2` | Non-confident |
| `3` | Classical circadian gene prefixes |

PNGs are written to `outdir/{fbase}plots_confident/` (or `_non_confident/`,
`_classic_circadian/`).
