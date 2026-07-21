# CHANGELOG — Unified scRNA-seq Pipeline

**Build date:** 2026-07-21
**Base:** `2026_nr4a1` (v10.0) — the most complete and balanced of the prior projects
**Headline change:** doublet detection is now **selectable**, with **DoubletFinder as the default**

---

## 1. Why this build exists

Four pipeline generations had drifted apart:

| Project | Role in this build |
|---|---|
| `github_pipeline` | The published/unified trail. Byte-identical to `2026_nr4a1` for scripts 00, 01, 02. |
| `2026_nr4a1` | Most complete and balanced. **Used as the architectural base.** Adds `06_annotation_unifier.R`, `08_cellchat.R`, and cleaned subannotation scripts. |
| `wu_project1` | Original DoubletFinder implementation. Source of the enhanced pK diagnostic plot. |
| `eithan_coffee` | Refined DoubletFinder implementation with constrained pK selection and rollback logging. **Source of the restored doublet module.** |

`github_pipeline` and `2026_nr4a1` were verified identical for `00_rlibs_installation.R`, `01_process_data.R`, and `02_global_annotation.R`; `2026_nr4a1` was the newer copy elsewhere. Nothing was lost by taking nr4a1 as the base.

---

## 2. File-by-file provenance

| File in `unified_pipeline/` | Origin | Changes made |
|---|---|---|
| `00_rlibs_installation.R` | `2026_nr4a1` v3.0 | → **v4.0.** Added `BPCells`; fixed `seuratWrappers` → `SeuratWrappers`; added hard dependency gate. |
| `01_process_data.R` | `2026_nr4a1` v10.0 | → **v11.0.** Doublet architecture rewritten (see §3). Everything else preserved. |
| `02_global_annotation.R` | `2026_nr4a1` v3.1 | Header tagged only. No logic changes. |
| `03_tcell_subannotation.R` | `2026_nr4a1` (`_clean` variant) | Renamed from `_clean`. Commented FeaturePlot updated to `DF_score`. |
| `04_colonocyte_subannotation.R` | `2026_nr4a1` (`_clean` variant) | Renamed from `_clean`. The older non-clean duplicate was dropped. |
| `05_macrophages_subannotation.R` | `2026_nr4a1` (`_clean` variant) | Renamed from `_clean`. Commented FeaturePlot updated to `DF_score`. |
| `06_annotation_unifier.R` | `2026_nr4a1` | Header tagged only. |
| `07_DE_and_two_way_ANOVA.R` | `2026_nr4a1` v1.1 | Header tagged only. Supersedes `github_pipeline/06_...`. |
| `08_cellchat.R` | `2026_nr4a1` | Header tagged only. |
| `copy_sc_data.sh` | `github_pipeline` + `eithan_coffee` | → **v2.0.** Merged and parameterized (see §5). |
| `cell_type_markers.csv` | `github_pipeline` | Unchanged. |
| `nr4a1_wholebody_ko_markers.txt` | `github_pipeline` | Unchanged. |
| `Nr4a1_s17_metadata.xlsx`, `Wu_metadata.xlsx` | `github_pipeline` | Unchanged — use as metadata templates. |
| `single-cell Standard Operating Procedure.docx` | `github_pipeline` | Unchanged. |

---

## 3. The doublet detection rewrite (Script 01)

### 3.1 What changed conceptually

v10.0 hardcoded `scDblFinder`. v11.0 introduces a **method switch** with **DoubletFinder restored as the default**, and normalizes both callers behind one interface so nothing downstream cares which ran.

```r
DOUBLET_METHOD         <- "DoubletFinder"   # "DoubletFinder" | "scDblFinder" | "both"
DOUBLET_CONSENSUS_RULE <- "intersect"       # only used when method == "both"
```

Both detectors are implemented as functions with an **identical contract** — `run_doubletfinder()` and `run_scdblfinder()` each return `$calls`, `$score`, and `$info`. The main loop dispatches and then writes a single standardized `Doublet_Status` column (`"Singlet"` / `"Doublet"`). **Scripts 02–08 therefore need no changes when you switch methods.**

### 3.2 The restored DoubletFinder module

Rebuilt from `eithan_coffee`, with additions:

- **paramSweep → summarizeSweep → find.pK** — the full sweep, as you had it.
- **Constrained pK selection** (from `eithan_coffee`): take the global BCmetric peak; if it falls outside `[DF_PK_RANGE_MIN, DF_PK_RANGE_MAX]` (default `[0.01, 0.15]`), fall back to the best peak *inside* the window; if no peak exists there, use `DF_PK_FALLBACK` (0.09). Every outcome is logged as `GLOBAL_PEAK_IN_RANGE` / `CONSTRAINED_TO_RANGE` / `HARDCODED_FALLBACK`.
- **Enhanced pK plot** (from `wu_project1`): grey BCmetric trace, red dashed line at the selection, red diamond on the selected point. **New:** a shaded band marking the plausible pK window, so a constrained selection is visually obvious.
- **Homotypic doublet adjustment** — *newly added, in neither original.* `modelHomotypic()` estimates the fraction of doublets that are same-cell-type and therefore undetectable, and shrinks `nExp` accordingly. Without this, DoubletFinder is forced to over-call heterotypic doublets to hit its `nExp` quota — a known cause of losing real transitional populations. Toggle: `DF_HOMOTYPIC_ADJUST` (default `TRUE`).
- **Small-sample guard** — *new.* Samples with too few cells for a stable PCA are skipped with all cells retained, instead of erroring out mid-run.

### 3.3 What was preserved from v10.0

- SCEVAN → DecontX dual-checkpoint system, untouched
- Rollback safeguard (`DOUBLET_ROLLBACK_THRESHOLD`), now applied to **both** methods
- Per-sample diagnostic UMAP of doublet scores
- scDblFinder path itself, fully intact and switchable

### 3.4 New: `"both"` mode

Runs both callers, then resolves per `DOUBLET_CONSENSUS_RULE`:

| Rule | Behavior |
|---|---|
| `intersect` | Doublet only if **both** agree — conservative, highest precision *(default)* |
| `union` | Doublet if **either** flags it — aggressive, highest recall |
| `DoubletFinder` | DF decides; scDblFinder logged only |
| `scDblFinder` | scDblFinder decides; DF logged only |

Also emits a **per-sample confusion matrix** plus **agreement % and Cohen's kappa**. Recommended workflow: run `"both"` once on a new dataset to see how far apart the callers are, then lock in a single method for production.

### 3.5 New: configuration validation

`validate_config()` runs before any sample is touched and hard-stops on impossible settings. Notably it catches `DOUBLET_ROLLBACK_THRESHOLD <= DOUBLET_RATE` — because DoubletFinder always calls approximately `DOUBLET_RATE` doublets by construction, that combination would silently trigger rollback on *every* sample and produce a run with no doublet removal at all.

### 3.5b New: automatic per-sample doublet rate

`DOUBLET_RATE <- "auto"` is now the default. Rather than applying one fixed rate to every sample, the rate is derived **per sample** from that sample's own recovered cell count, using the 10x linear multiplet model:

```
multiplet rate ≈ 0.8% per 1,000 cells recovered
```

| Cells recovered | Auto rate |
|---|---|
| 5,000 | 4.0% |
| 10,000 | 8.0% |
| 16,000 | 12.8% |
| 20,000 | 16.0% |

This replaces the inherited hardcoded `0.16`, which was only correct for samples recovering ~20,000 cells. Applied to a 6,000-cell sample it implied a 16% doublet quota against a true expectation near 4.8% — under DoubletFinder that is a hard quota, so roughly 11% of that sample's real cells were being deleted to fill it.

Controls:

| Parameter | Default | Purpose |
|---|---|---|
| `DOUBLET_RATE` | `"auto"` | `"auto"` or a fixed number |
| `DOUBLET_RATE_PER_1K` | `0.008` | The 10x slope |
| `DOUBLET_RATE_AUTO_MIN` | `0.008` | Floor for tiny samples |
| `DOUBLET_RATE_AUTO_MAX` | `0.20` | Cap; also keeps auto rates below the rollback threshold |

The rate used and its provenance (`AUTO` / `AUTO_CAPPED` / `AUTO_FLOORED` / `FIXED`) are recorded per sample in the doublet log, and the console prints the observed range across samples at the end of the step. A capped sample raises a warning, since that usually indicates an overloaded lane worth checking.

### 3.6 New metadata columns

| Column | Written when |
|---|---|
| `Doublet_Status` | Always — the standardized final call |
| `DF_score`, `DF_class` | DoubletFinder ran (pANN score + raw call) |
| `scDblFinder_score`, `scDblFinder_class` | scDblFinder ran |

---

## 3A. New scripts: 09 and 10

### `09_cell_scores.R` — potency, stemness, entropy

Three independent estimates of differentiation state, so no conclusion rests on a single method:

| Method | Type | Output | Scale |
|---|---|---|---|
| **CytoTRACE 2** | Supervised deep learning | `CytoTRACE2_Score`, `CytoTRACE2_Potency` (6 categories), `CytoTRACE2_Relative`, plus pre-KNN variants | **Absolute**, comparable across datasets |
| **CytoTRACE v1** | Unsupervised (gene counts) | `CytoTRACE1_Score`, `_Rank`, `_GCS` | Relative to this run only |
| **Entropy** | Native, no dependencies | `entropy_shannon`, `entropy_normalized`, `gene_counts_score` | Relative |

Key design points:

- **Graceful degradation.** Each method is checked for availability up front; a missing package is skipped with a pointer to `INSTALL_NOTES.md` and the rest still run. The entropy metrics have zero external dependencies, so you always get a potency readout.
- **Raw counts enforced.** CytoTRACE 2 must not receive log-transformed data; the script pulls from `counts`, never `data`.
- **Per-sample by default** (`CT2_RUN_PER_SAMPLE <- TRUE`), following the CytoTRACE 2 FAQ — the KNN post-processing borrows across neighbours and is distorted by batch effects.
- **Entropy is depth-normalized** (`H / log(n_detected)`). Raw Shannon entropy is dominated by sequencing depth; without normalization you would mostly be measuring library size.
- **Concordance check.** All methods should correlate *positively* (higher = less differentiated). The script computes a Spearman matrix, writes a heatmap, and raises a warning on any negative pair — that is a red flag, not a curiosity.
- **Honest statistics.** Group comparisons are non-parametric with BH correction, and every output carries an explicit caveat that cell-level tests treat cells as independent replicates when they are not. Confirm key findings at the sample level.

Outputs: `<PROJECT>_with_cell_scores.rds`, `cell_scores/` plots, `cell_scores_summary.xlsx`, and a gzipped per-cell table.

### `10_trajectory_cellrank.R` — directed trajectories

Script 09 asks how potent each cell is; Script 10 asks where it is going. Uses CellRank's **CytoTRACEKernel**, which derives directionality from the CytoTRACE gene-count signal — so it works **without RNA velocity**, which matters because standard Cell Ranger output has no spliced/unspliced layers.

- Runs Python (scanpy + CellRank) through `reticulate`; the environment is created once by hand, never auto-built from inside the script.
- **Verifies every Python module before doing any work**, so misconfiguration fails in seconds rather than after loading a large object.
- Seurat → AnnData conversion is done manually rather than via SeuratDisk, which is unmaintained and breaks on Seurat v5.
- Reuses the **Harmony embedding** by default, so CellRank does not mistake sample-of-origin for a differentiation axis.
- `SUBSET_TO` restricts to one lineage — with a loud warning if left `NULL`, because transitions modelled between unrelated lineages (T cells ↔ colonocytes) are meaningless.

Outputs: pseudotime (`cellrank_ct_pseudotime`), fate probabilities per terminal state, `fate_dominant`/`fate_confidence`, plots, and an optional `.h5ad` for Python follow-up.

### `INSTALL_NOTES.md`

Install guide for the packages that actually break, with per-package failure tables: CytoTRACE 2 (the `subdir` trap), CytoTRACE v1 (stale deps — skippable), DoubletFinder (renamed functions), BPCells, SCEVAN, CellChat, monocle3, and the full reticulate/conda setup for CellRank including the "restart R, reticulate binds once" gotcha. Ends with a copy-paste verification block and a required-vs-skippable table.

### 3.7 New outputs

All diagnostics now land in one place: **`seurat_output/doublet_diagnostics/`**

- `<sample>_DF_pk_plot.png` — pK selection with plausible-range band
- `<sample>_scDbl_score_distribution.png` — score histogram
- `<sample>_doublet_umap.png` — doublet score on UMAP
- `<sample>_doublet_concordance.png` — confusion matrix (`"both"` mode only)
- `seurat_output/doublet_detection_log_<METHOD>.xlsx` — master audit table: method, pK used vs. global peak, nExp raw vs. homotypic-adjusted, doublet fraction, agreement/kappa, and rollback action per sample

---

## 4. Script 00 fixes

1. **`BPCells` added.** Script 01 calls `library(BPCells)` but it appeared in no install list — `USE_BPCELLS <- TRUE` would have failed at load.
2. **`seuratWrappers` → `SeuratWrappers`.** The namespace is capitalized, so `requireNamespace("seuratWrappers")` never matched an installed package and forced a GitHub rebuild on every single run.
3. **Hard dependency gate.** The script now `stop()`s with an explicit list if a Script 01 dependency is missing, rather than printing a status table that is easy to skim past.

---

## 5. `copy_sc_data.sh` v2.0

Merged from both originals and parameterized.

- **Bug fixed:** `eithan_coffee`'s version had `path2` commented out but still ran `cp -rv ${path2}` on the empty variable — a cp error on every sample. Optional files are now behind explicit toggles (`COPY_PROBE`, `COPY_SUMMARY`).
- Configurable `SOURCE_ROOT` / `POOLS` / `DEST_DIR`; preset pool lists for both studies are in the comments.
- Iterates sample directories directly instead of parsing `ls`, so unusual sample names are safe.
- Existence checks, `OVERWRITE=no` re-run safety, `DRY_RUN=1` preview mode.
- **Exits non-zero if any required H5 is missing**, so a silently-absent sample cannot slip into the R pipeline.

---

## 6. Migration notes

- **Object compatibility:** objects produced by the old v10.0 script still work downstream — `Doublet_Status` has the same name and meaning.
- **Re-running:** the checkpoint system is unchanged, so existing `_scevan_processed.rds` / `_decontx_dblt_processed.rds` files are still honored. **To re-run doublet detection with the new method, delete the Checkpoint 2 files** (`*_decontx_dblt_processed.rds`); Checkpoint 1 is retained so SCEVAN will not re-run.
- **Paths:** `ROOT_PATH` in Script 01 still points at `/home/ssromerogon/2026_nr4a1_ack/r_process`. Update it, along with `PROJECT_NAME` and `METADATA_FILE`, for a new study.
- **Probe settings:** `ADD_PROBE_DATA <- TRUE` with Nr4a1 exon probes is carried over from the nr4a1 study. Set `FALSE` for studies without a probe assay.

---

## 7. Not carried over

- `eithan_coffee/03_de_and_two_way_anova_eithan.R` — you chose the full nr4a1 set; `07_DE_and_two_way_ANOVA.R` (v1.1, SplineDV fix) supersedes it.
- `eithan_coffee/02_annotate_data_SG.R` — superseded by `02_global_annotation.R` (CSV-driven, weighted pre-scoring).
- `wu_project1/02_cell_type_annotations.R` — project-specific hardcoded annotations.
- `2026_nr4a1/04_colonocyte_subannotation.R` (non-clean) — the `_clean` variant was kept.
