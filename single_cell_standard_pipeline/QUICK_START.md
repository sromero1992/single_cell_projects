# Single-Cell RNA-seq Pipeline — Quick Start (Unified Build)

> Full methodology and parameter rationale: **single-cell Standard Operating Procedure.docx** in this folder.
> What changed relative to the older project folders: **CHANGELOG.md**.

**This build's headline feature:** doublet detection is selectable via one variable, with **DoubletFinder restored as the default**.

---

## Prerequisites

1. R ≥ 4.3 with RStudio (or compatible IDE).
2. Run `00_rlibs_installation.R` **once per machine**. It now hard-stops if anything Script 01 needs is missing, so a clean exit means you are genuinely ready.
3. A filled-in metadata `.xlsx` (format below). `Nr4a1_s17_metadata.xlsx` and `Wu_metadata.xlsx` are usable templates.
4. H5 files staged by `copy_sc_data.sh`.
5. For scripts 09/10 only: see **INSTALL_NOTES.md** — CytoTRACE 2 needs a `subdir` argument, and Script 10 needs a one-time Python environment.

---

## Run Order

| Step | Script | What it does | Typical runtime |
|---|---|---|---|
| — | `copy_sc_data.sh` | Stage Cell Ranger H5s into `h5_files/` | 5–30 min |
| 0 | `00_rlibs_installation.R` | Install all dependencies | 15–60 min (first run only) |
| 1 | `01_process_data.R` | Load H5 → SCEVAN → DecontX → **doublets** → QC → Harmony | 1–4 h per 10 samples |
| 2 | `02_global_annotation.R` | Weighted pre-scoring → broad cell type annotation | Interactive, ~30 min |
| 3 | `03_tcell_subannotation.R` | T cell sub-clustering and annotation | Interactive |
| 4 | `04_colonocyte_subannotation.R` | Colonocyte sub-clustering and annotation | Interactive |
| 5 | `05_macrophages_subannotation.R` | Macrophage sub-clustering and annotation | Interactive |
| 6 | `06_annotation_unifier.R` | Merge all sub-annotations into one final object | Fast |
| 7 | `07_DE_and_two_way_ANOVA.R` | Wilcoxon DE + two-way ANOVA + SplineDV + GSEA + Enrichr | 30 min – 2 h |
| 8 | `08_cellchat.R` | Cell–cell communication analysis | 30 min – 2 h |
| 9 | `09_cell_scores.R` | CytoTRACE 2 + CytoTRACE v1 + entropy potency scoring | 30 min – 2 h |
| 10 | `10_trajectory_cellrank.R` | CellRank CytoTRACEKernel trajectories (needs Python) | 30 min – 1 h |

Scripts run **in order** — each reads the `.rds` output of the previous. Scripts 03–05 are independent of each other and can run in any order, but all must finish before `06`.

---

## Staging data

Edit the CONFIGURATION block of `copy_sc_data.sh` (`SOURCE_ROOT`, `POOLS`, toggles), then:

```bash
chmod +x copy_sc_data.sh
DRY_RUN=1 ./copy_sc_data.sh    # preview
./copy_sc_data.sh              # copy for real
```

It exits non-zero and lists the offenders if any required H5 is missing. Every `SampleID` in your metadata must have an exactly-matching folder in `h5_files/`.

---

## Doublet Detection — the main new control

All in **Script 01, Section 1.4**.

```r
DOUBLET_METHOD <- "DoubletFinder"   # "DoubletFinder" | "scDblFinder" | "both"
```

| Value | When to use |
|---|---|
| `"DoubletFinder"` | **Default.** paramSweep + constrained pK. Conservative; preserves transitional and proliferating populations. |
| `"scDblFinder"` | Faster, fully automatic, no pK. Better in published benchmarks but was over-aggressive on this lab's data. |
| `"both"` | Run once on a new dataset to see how far apart the two callers are, then commit to one. |

### Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `DOUBLET_RATE` | `"auto"` | **Computed per sample** from that sample's recovered cell count (~0.8% per 1,000 cells). Set a number to override. Under DoubletFinder this is a **hard quota**, not a ceiling. |
| `DOUBLET_RATE_AUTO_MAX` | `0.20` | Cap on the automatic estimate; also keeps it below the rollback threshold. |
| `DOUBLET_ROLLBACK_THRESHOLD` | `0.25` | Must stay **above** `DOUBLET_RATE` or rollback fires on every sample. Validated at startup. |
| `DF_PK_RANGE_MIN` / `MAX` | `0.01` / `0.15` | Plausible pK window. If the BCmetric global peak lands outside it, the best in-window peak is used instead. |
| `DF_PK_FALLBACK` | `0.09` | Used only if no peak exists in the window at all. |
| `DF_HOMOTYPIC_ADJUST` | `TRUE` | Shrinks `nExp` by the estimated same-cell-type (undetectable) doublet fraction. Leave on. |
| `DOUBLET_CONSENSUS_RULE` | `"intersect"` | `"both"` mode only. `intersect` = conservative, `union` = aggressive. |

### What to check after Script 01

Everything lands in `seurat_output/doublet_diagnostics/` plus the master log `seurat_output/doublet_detection_log_<METHOD>.xlsx`.

1. **`Action` column** — any `ROLLBACK_APPLIED` or `ERROR_SKIPPED` rows mean those samples kept *all* cells.
2. **`Note` column** — `HARDCODED_FALLBACK` means the sweep found no usable peak; inspect that sample's pK plot.
3. **`pK_Used` vs `pK_Global_Peak`** — if they differ, constraint logic intervened. Expected occasionally; suspicious if it happens on every sample.
4. **`Target_Rate_Pct` and `Rate_Source`** — the doublet rate actually applied per sample. `AUTO_CAPPED` means the lane looked overloaded and the estimate was clamped; worth checking.
5. **`<sample>_DF_pk_plot.png`** — the selected pK should sit on a real peak inside the shaded band, not on a noisy shoulder.

> **Re-running doublets with a different method:** delete the Checkpoint 2 files (`*_decontx_dblt_processed.rds`) inside `seurat_output/scevan_per_sample_results/<SampleID>/`. Checkpoint 1 is kept, so SCEVAN will not re-run.

---

## Parameters You Touch Every New Experiment

### Script 01, Section 1.1

| Parameter | What to change |
|---|---|
| `PROJECT_NAME` | Short alphanumeric tag for this study |
| `ROOT_PATH` | Top-level project directory (currently `/home/ssromerogon/2026_nr4a1_ack/r_process`) |
| `METADATA_FILE` | Path to your sample metadata `.xlsx` |
| `H5_DIR` | Should point at the `h5_files/` produced by `copy_sc_data.sh` |

### Script 01, QC and toggles

| Parameter | Default | Notes |
|---|---|---|
| `POST_MIN_GENES` | 500 | Raise if low-quality cells persist |
| `POST_MAX_MT` | 5.0 | Adjust after viewing QC violins; standard range 5–20 |
| `POST_MAX_UMIS` | 100000 | Tune from QC violins |
| `RUN_SCEVAN` | TRUE | Enable if **any** sample is tumor — run on all samples, normals serve as CNA baseline. FALSE if the whole study is normal tissue. |
| `RUN_DECONTX` | TRUE | Recommended |
| `USE_BPCELLS` | FALSE | Set TRUE for large cohorts (>8 samples) to limit RAM. Requires the BPCells package. |
| `ADD_PROBE_DATA` | TRUE | **Set FALSE** unless this run used a probe assay. Currently configured for Nr4a1 exon probes. |

### Scripts 07 & 08

| Parameter | What to change |
|---|---|
| `CONDITION_COLUMN` | Metadata column defining the comparison (e.g. `"Condition"`) |
| `GROUP_1` / `GROUP_2` | `GROUP_1` = KO/treatment, `GROUP_2` = WT/control — positive log2FC means up in `GROUP_1` |
| `CONDITION_LEVELS` | Plot axis order: **control first**, e.g. `c("WT", "KO")` |

> **Convention:** `GROUP_1 = KO` so positive log2FC always means *up in KO*; `CONDITION_LEVELS = c("WT","KO")` so WT is always the visual baseline. These are intentionally independent.

---

## Cell Scores — Script 09

Three independent readouts of differentiation state, so no conclusion depends on one method.

| Method | Toggle | Output columns | Scale |
|---|---|---|---|
| CytoTRACE 2 | `RUN_CYTOTRACE2` | `CytoTRACE2_Score`, `CytoTRACE2_Potency`, `CytoTRACE2_Relative` | **Absolute** — comparable across datasets |
| CytoTRACE v1 | `RUN_CYTOTRACE1` | `CytoTRACE1_Score`, `_Rank`, `_GCS` | Relative to this run |
| Entropy | `RUN_ENTROPY` | `entropy_shannon`, `entropy_normalized`, `gene_counts_score` | Relative |

**All three point the same way: higher = less differentiated.**

Key parameters:

| Parameter | Default | Notes |
|---|---|---|
| `CT2_SPECIES` | `"mouse"` | `"human"` triggers internal orthology mapping |
| `CT2_RUN_PER_SAMPLE` | `TRUE` | Recommended — the KNN smoothing step is distorted by batch effects |
| `CT2_SLOT` | `"counts"` | **Must be raw/CPM counts.** Never log-transformed |
| `CT2_NCORES` | `4` | Drop to 1–2 on machines under 16 GB RAM |
| `CT2_USE_PREKNN` | `FALSE` | Set `TRUE` if rare populations (≤5 cells) are the point of the study |
| `ENTROPY_NORMALIZE` | `TRUE` | Divides by `log(n_detected)`. Leave on — raw entropy mostly measures sequencing depth |

**What to check:** open `cell_scores/method_concordance_heatmap.png` first. Every pairwise Spearman correlation should be **positive**. A negative pair means the methods disagree on direction — investigate before interpreting anything.

The entropy metrics have no external dependencies, so Script 09 produces a usable result even if both CytoTRACE packages fail to install.

> **On the statistics:** group comparisons in `cell_scores_summary.xlsx` treat individual cells as independent replicates. They are not — cells from one animal are correlated, so p-values are anti-conservative. Use them to rank effects, then confirm anything important with a sample-level test (per-sample medians, n = animals).

---

## Trajectories — Script 10

Uses CellRank's **CytoTRACEKernel**, which derives direction from the gene-count signal and therefore needs **no RNA velocity** — important, since standard Cell Ranger output has no spliced/unspliced layers.

**Requires a Python environment.** One-time setup in INSTALL_NOTES.md section 5:

```bash
conda create -n cellrank_env python=3.11 -y
conda activate cellrank_env
pip install "cellrank>=2.0" scanpy scvelo anndata igraph leidenalg
```

| Parameter | Default | Notes |
|---|---|---|
| `PYTHON_ENV` | `"cellrank_env"` | Must match the env you created |
| `SUBSET_TO` | `NULL` | **Set this.** Restrict to one lineage — see below |
| `USE_HARMONY_EMBEDDING` | `TRUE` | Stops CellRank reading batch structure as a differentiation axis |
| `CT_IMPUTE_METHOD` | `"moments"` | Smoothing for the kernel; `"cellmapper"` or `"none"` also available |
| `MAX_CELLS` | `30000` | Downsample above this |
| `N_STATES` / `N_TERMINAL_STATES` | `5` / `3` | Lower these if state estimation fails |

> ⚠️ **Set `SUBSET_TO` to a single differentiation hierarchy.** Trajectory inference across unrelated lineages models transitions between cell types that cannot interconvert — T cells "becoming" colonocytes — and the output is meaningless. The script warns loudly if you leave it `NULL`.

Outputs `cellrank_ct_pseudotime` (low = less differentiated), per-terminal-state fate probabilities, and `fate_dominant`/`fate_confidence`.

---

## Metadata File Format

**Required:** `SampleID` — must exactly match the subfolder name inside `H5_DIR`.

```
SampleID | Condition | Genotype | Diet | Sex | Timepoint | BatchID
S01_WT_M |    WT     |    WT    | ctrl |  M  |   Week4   |  Batch1
S02_KO_M |    KO     | Nr4a1_KO | ctrl |  M  |   Week4   |  Batch1
```

Every column is attached to every cell from that sample. Auto-created: `group` = `Condition_Sex` (if both exist); any `ColA_ColB` pattern is built on demand when referenced downstream.

---

## Marker Gene CSV Format

`cell_type_markers.csv` — one row per cell type/marker combination:

| Column | Description |
|---|---|
| `cell_type` | Display name (e.g. `T cells`) |
| `parent_cell_type` | Broad type this sub-type belongs to (same as `cell_type` for broad markers) |
| `tier` | `broad` or `sub` |
| `markers` | Pipe-separated symbols: `Cd3e\|Cd3d\|Cd3g` |
| `sub_resolution` | Leiden resolution for sub-clustering |

---

## Output Folder Structure

```
seurat_output/
├── <PROJECT>_processed_for_annotation.rds    # Script 01 → 02 handoff
├── <PROJECT>_final_annotated.rds             # Script 06 → 07/08 handoff
├── doublet_detection_log_<METHOD>.xlsx       # doublet audit trail
├── doublet_diagnostics/                      # pK plots, score histograms,
│   │                                         #   doublet UMAPs, concordance
│   ├── <sample>_DF_pk_plot.png
│   ├── <sample>_doublet_umap.png
│   └── <sample>_doublet_concordance.png      # "both" mode only
├── scevan_per_sample_results/                # per-sample SCEVAN + checkpoints
│   └── <SampleID>/
│       ├── <SampleID>_scevan_processed.rds        # Checkpoint 1
│       └── <SampleID>_decontx_dblt_processed.rds  # Checkpoint 2
├── QC_plots/                                 # pre/post-filter violins
├── annotation_plots/                         # UMAPs, dotplots, heatmaps
├── differential_expression/
├── cellchat/
├── <PROJECT>_with_cell_scores.rds            # Script 09 → 10 handoff
├── cell_scores/                              # Script 09
│   ├── method_concordance_heatmap.png        #   CHECK THIS FIRST
│   ├── umap_<score>.png
│   ├── box_<score>_by_celltype.png
│   ├── composition_potency_by_<group>.png
│   ├── cell_scores_summary.xlsx
│   └── cell_scores_per_cell.csv.gz
├── <PROJECT>_with_trajectory.rds             # Script 10 output
└── trajectory/                               # Script 10
    ├── umap_cellrank_ct_pseudotime.png
    ├── umap_fate_dominant.png
    ├── pseudotime_by_celltype.png
    ├── trajectory_summary.xlsx
    └── adata_cellrank.h5ad                   # for Python follow-up
```

---

## Common Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Script 00 stops with "MISSING CRITICAL PACKAGES" | A dependency failed to build | Install the listed packages manually; for `DoubletFinder` confirm `remotes` and a compiler toolchain are present |
| Script 01 stops immediately with "CONFIGURATION ERRORS" | Invalid Section 1.4 combination | Read the message — most often `DOUBLET_ROLLBACK_THRESHOLD <= DOUBLET_RATE` |
| Every sample shows `ROLLBACK_APPLIED` | `DOUBLET_RATE` too close to the rollback threshold | Use `DOUBLET_RATE <- "auto"`, or lower it / raise the threshold |
| Script 09 skips CytoTRACE 2 | Package not installed | See INSTALL_NOTES.md section 1 — the `subdir` argument is required |
| Script 09 concordance heatmap shows negative rho | Methods disagree on direction | Investigate before interpreting; check the counts layer is not log-transformed |
| Script 10 binds the wrong Python | reticulate binds once per session | **Restart R**, then set `PYTHON_ENV` before anything imports Python |
| Script 10: `CytoTRACEKernel` not found | CellRank v1 installed | `pip install --upgrade "cellrank>=2.0"` |
| Doublet removal deletes real biology | Method too aggressive | Switch to `"DoubletFinder"`, confirm `DF_HOMOTYPIC_ADJUST <- TRUE`, lower `DOUBLET_RATE` |
| `HARDCODED_FALLBACK` on many samples | BCmetric curve has no in-range peak | Inspect pK plots; widen `DF_PK_RANGE_MAX` if the real peak sits just outside |
| paramSweep extremely slow | Expected — it is the slow step | Use `"scDblFinder"` for a fast first pass, or run overnight |
| OOM on >10 samples | RAM exhausted during SCEVAN | Restart — checkpoints resume automatically; consider `USE_BPCELLS <- TRUE` |
| `library(BPCells)` fails | BPCells not installed | Set `USE_BPCELLS <- FALSE`, or install via Script 00 |
| Enrichr returns nothing | No internet or wrong DB name | Verify with `enrichR::listEnrichrDbs()` |

---

*Full parameter reference: single-cell Standard Operating Procedure.docx*
