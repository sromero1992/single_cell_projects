# TamuScDSC — modular single-cell preprocessing, QC & doublet detection

**Status:** design + scaffold (v0.1). Migrated from `01_process_data_new.R` and
`01b_cluster_doublet_review.R`.

**Goal.** Turn the preprocessing / QC / doublet stage from a single procedural
script into an installable R package whose steps are *atomic, order-free, and
re-runnable*, and which accepts three different starting points without forking
the code. Scripts 02–10 stay as scripts; they just `library(TamuScDSC)` for the
shared markers / metadata / utilities.

---

## 1. The one idea: normalize every input to a *sample list*

The single design decision that makes everything else fall out is to convert all
three input scenarios into one internal representation at the boundary — a
**named list of per-sample Seurat objects** — and have every step operate on that.

```
Scenario 1: one big merged Seurat object   --SplitObject(split.by=SampleID)-->  list
Scenario 2: per-sample .rds checkpoints     --readRDS each-------------------->  list
Scenario 3: preloaded list of Seurat objs   --already a list------------------>  list
```

One generic does the conversion:

```r
samples <- as_sample_list(x, sample_col = "SampleID")   # x = paths | Seurat | list
```

and `merge_samples(samples)` collapses the list back into one object when you are
ready to integrate. "Where the data came from" is now completely decoupled from
"what you do to it." All three scenarios are first-class — the generic dispatches
on the class of `x`.

## 2. Atomic steps that read/write a fixed metadata schema

Each stage is a standalone function: object(s) in, object(s) out. A step never
assumes what ran before it and never deletes cells — it only *writes agreed
metadata columns*. The schema (see `R/schema.R`) is the contract between steps:

| Column                | Written by                | Meaning                                  |
|-----------------------|---------------------------|------------------------------------------|
| `percent_mt`          | `qc_light` / `qc_stringent` | mitochondrial %                        |
| `qc_pass_light`       | `qc_light`                | passes loose per-sample QC               |
| `qc_pass_stringent`   | `qc_stringent`            | passes stringent merged QC               |
| `scevan_class`        | `run_scevan`             | tumor / normal / filtered                |
| `DF_score`,`DF_class` | `detect_doublets`        | DoubletFinder pANN + call                |
| `scDblFinder_score`,`scDblFinder_class` | `detect_doublets` | scDblFinder prob + call        |
| `doublet_consensus`   | `detect_doublets`        | Singlet/Doublet after the consensus rule |
| `cluster_dbl_flag`    | `cluster_doublet_review` | keep / flagged_doublet_cluster           |

Because detection only *labels*, a single function — `apply_doublet_action()` —
is the **only** thing that ever removes cells. That is what turns "remove both /
remove consensus / keep and flag / let 01b decide then remove" into an *argument*
rather than a fork in the code.

## 3. Ordering freedom: you compose, the package does not decide

The you-vs-colleague disagreement (doublets per-sample-before-merge vs.
doublets-after-stringent-merged-QC) is resolved by not baking an order in. The
package ships **both** as named recipes, and each is just a short pipe over the
same atoms:

```r
## Recipe A — per-sample, doublets before merge (your current approach)
data <- recipe_per_sample(input, metadata)
# == as_sample_list(input) |> attach_metadata(metadata) |>
#    qc_light() |> run_decontx() |> detect_doublets(method="both") |>
#    merge_samples() |> qc_stringent() |> integrate_data()

## Recipe B — stringent merged QC first, doublets after (colleague's approach)
data <- recipe_merged_stringent(input, metadata)
# == as_sample_list(input) |> attach_metadata(metadata) |> merge_samples() |>
#    qc_stringent() |> detect_doublets(method="both", by="SampleID") |>
#    cluster_doublet_review(action="flag") |> integrate_data()
```

`detect_doublets()` is level-agnostic: given a **list** it maps per sample; given
a **merged** object it splits by `by=` internally (or runs global if `by=NULL`).
That is what lets the same step run before *or* after the merge.

## 4. Re-running one sample and grafting it back

Every join is keyed on cell barcode, never on row position. Two helpers give you
the "extract a sample, redo a step, put it back" workflow:

```r
sub  <- subset_samples(data, ids = "WT_3")
sub  <- detect_doublets(sub, method = "both", rate = 0.08)          # new params
data <- graft_meta(data, sub, cols = c("DF_score","scDblFinder_score",
                                       "doublet_consensus"))
```

`graft_meta()` writes the recomputed columns back into the parent by barcode, so
you iterate on a problem sample without touching the other eleven.

## 5. Markers as package data, hard-coded but extensible

The default marker set ships inside the package and is edited through an API, not
by touching code:

```r
m <- get_markers()                                   # built-in named list
m <- add_markers(m, "Tuft", c("Dclk1","Trpm5","Gfi1b"))
m <- load_markers("my_project_markers.csv")          # override / extend from file
```

Scripts 02–05 call `get_markers()` instead of each carrying their own copy — one
source of truth, per-project extensions layered on top.

## 6. Provenance without losing freedom

`run_pipeline()` (and each step) stamps `obj@misc$TamuScDSC_provenance` with the step
name, timestamp, and parameters. `provenance(obj)` prints what happened to an
object. This is how total ordering freedom stays auditable.

---

## 7. Package layout

```
TamuScDSC/
  DESCRIPTION            # deps, metadata
  NAMESPACE              # exports (regenerate with devtools::document())
  R/
    schema.R             # column-name contract + validators
    utils.R              # %||%, provenance stamp, small helpers
    ingest.R             # as_sample_list() S3, read_checkpoints(), merge_samples()
    metadata.R           # clean_header(), norm_key(), resolve_sample_col(), attach_metadata()
    qc.R                 # qc_light(), qc_stringent()
    doublets.R           # resolve_doublet_rate(), run_doubletfinder(), run_scdblfinder(),
                         #   report_concordance(), detect_doublets(), apply_doublet_action()
    scevan.R             # run_scevan()
    decontx.R            # run_decontx()
    cluster_review.R     # cluster_doublet_review()  (the 01b logic)
    integrate.R          # integrate_data()  (normalize + HVG + PCA + Harmony/IntegrateLayers + UMAP)
    markers.R            # get_markers(), add_markers(), load_markers()
    subset.R             # subset_samples(), graft_meta()
    recipes.R            # recipe_per_sample(), recipe_merged_stringent()
    pipeline.R           # run_pipeline(), provenance()
  inst/extdata/
    cell_type_markers.csv
  vignettes/
    three-scenarios.Rmd  # one worked example per input type
  tests/testthat/        # locks the schema + barcode-join behaviour
```

## 8. How the existing scripts map in

| Old location (`01_process_data_new.R`)        | New home                          |
|-----------------------------------------------|-----------------------------------|
| `validate_config()`                           | folded into each step's arg checks + `doublet_params()` |
| `resolve_doublet_rate()`                      | `R/doublets.R` (params as args)   |
| `run_doubletfinder()` / `run_scdblfinder()`   | `R/doublets.R` (verbatim logic)   |
| `report_concordance()`                        | `R/doublets.R`                    |
| `clean_header()` / `norm_key()`               | `R/metadata.R`                    |
| STEP 2.1a per-sample loop                     | `qc_light` + `run_decontx` + `detect_doublets` mapped over the list |
| STEP 2.2 merge                                | `merge_samples()`                 |
| STEP 2.3 post-merge QC                        | `qc_stringent()`                  |
| STEP 2.4 normalize/reduce/integrate           | `integrate_data()`                |
| `01b_cluster_doublet_review.R`                | `cluster_doublet_review()`        |

Global `CONFIG` constants become **function arguments with defaults**, grouped
into small param-list constructors (`doublet_params()`, `qc_params()`) so a
project still has one place to set everything, but steps are callable in isolation.

## 9. What is intentionally NOT in the package

Annotation (02–06), DE/DV (07), CellChat (08), cell scores (09), trajectories
(10) stay as standalone scripts. They gain a `library(TamuScDSC)` line for markers /
metadata / utilities but keep their per-project freedom. The package deliberately
owns only the stage that benefits from being locked down: ingestion, QC, doublets,
integration, and the shared marker/metadata plumbing.
