# TamuScDSC

Modular, installable R package for the **preprocessing / QC / doublet** stage of
the lab single-cell pipeline. Migrated from `01_process_data_new.R` and
`01b_cluster_doublet_review.R`. See [`ARCHITECTURE.md`](ARCHITECTURE.md) for the
full design.

## Why

The preprocessing stage is the part that needs the *least* freedom and the most
consistency across datasets. This package makes each step atomic, order-free, and
re-runnable, and accepts three starting points interchangeably.

## Install

```r
# from the package directory
devtools::install("TamuScDSC")          # or: R CMD INSTALL TamuScDSC
# regenerate NAMESPACE + man/ from roxygen if you edit exports:
devtools::document("TamuScDSC")
```

Heavy dependencies (DoubletFinder, scDblFinder, harmony, celda, SCEVAN) are
*Suggests*: the package installs without them, and each step errors with a clear
install hint only if you actually call it.

## 30-second tour

```r
library(TamuScDSC)

# 1) any input -> a named per-sample list
samples <- as_sample_list(x)              # x = paths | Seurat object | list

# 2) run a shipped recipe ...
data <- recipe_per_sample(samples, metadata = meta_df)        # doublets pre-merge
data <- recipe_merged_stringent(samples, metadata = meta_df)  # doublets post-QC

# 3) ... or compose the atoms yourself
data <- samples |>
  qc_light() |> run_decontx() |> detect_doublets(method = "both") |>
  merge_samples() |> qc_stringent() |> integrate_data(method = "RunHarmony")

# 4) review clusters, remove, re-run a sample, graft back
data <- cluster_doublet_review(data, action = "flag", out_dir = "review")
data <- apply_doublet_action(data, action = "remove", target = "cluster_flag")
sub  <- detect_doublets(subset_samples(data, "WT_3"), rate = 0.08)
data <- graft_meta(data, sub, cols = c("DF_score","scDblFinder_score","doublet_consensus"))
```

## Key functions

| Area        | Functions |
|-------------|-----------|
| Ingest      | `as_sample_list()`, `read_checkpoints()`, `merge_samples()` |
| Metadata    | `attach_metadata()`, `clean_header()`, `norm_key()`, `resolve_sample_col()` |
| QC          | `qc_light()`, `qc_stringent()`, `qc_params()` |
| Doublets    | `detect_doublets()`, `apply_doublet_action()`, `run_doubletfinder()`, `run_scdblfinder()`, `report_concordance()`, `resolve_doublet_rate()`, `doublet_params()` |
| Review      | `cluster_doublet_review()`, `cluster_review_params()` |
| Other steps | `run_decontx()`, `run_scevan()`, `integrate_data()` |
| Markers     | `get_markers()`, `add_markers()`, `load_markers()`, `markers_as_list()` |
| Subset      | `subset_samples()`, `graft_meta()` |
| Orchestrate | `recipe_per_sample()`, `recipe_merged_stringent()`, `run_pipeline()`, `provenance()` |

## Contract

All steps read/write a fixed metadata schema (`TamuScDSC_schema()`), and only
`apply_doublet_action()` ever removes cells. Every join is barcode-keyed. Scripts
02–10 stay as standalone scripts and just `library(TamuScDSC)` for markers, metadata
and utilities.
