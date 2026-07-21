# Single-Cell RNA-seq Pipeline — Quick Start Guide

> For full methodology, parameter rationale, and output file descriptions see the
> **Single-cell Standard Operating Procedure.docx** in this folder.

---

## Prerequisites

1. R ≥ 4.3 with RStudio (or a compatible IDE).
2. All packages installed — run `00_rlibs_installation.R` **once per machine**
   before doing anything else.
3. A filled-in metadata Excel file (see Section 3 of the SOP or the template in
   Script 01 for the required columns).

---

## Run Order

| Step | Script | What it does | Typical runtime |
|------|--------|-------------|-----------------|
| 0 | `00_rlibs_installation.R` | Install all dependencies | 15–60 min (first run only) |
| 1 | `01_process_data.R` | Load H5 → QC → doublet removal → Harmony | 1–4 h per 10 samples |
| 2 | `02_annotate_data.R` | Weighted pre-scoring → manual annotation → sub-clustering | Interactive, ~30 min |
| 3 | `03_differential_expression.R` | Wilcoxon DE + SplineDV + GSEA + Enrichr ORA | 30 min – 2 h |
| 4 | `04_cell_scoring.R` | AUCell pathway scoring + pseudo-bulk t-tests | 20–60 min |

Scripts **must be run in order**. Each script reads the `.rds` output of the previous one.

---

## The 10 Parameters You Touch Every New Experiment

Edit these in **Part 1** of each script before running. Everything else can stay at its default.

### Scripts 01–04 (shared)

| Parameter | Where | What to change |
|-----------|-------|----------------|
| `PROJECT_NAME` | All scripts, Section 1.1 | Short alphanumeric tag for this study |
| `ROOT_PATH` | All scripts, Section 1.1 | Top-level directory for this project |
| `METADATA_FILE` | Script 01, Section 1.1 | Path to your sample metadata `.xlsx` |

### Script 01

| Parameter | Default | Notes |
|-----------|---------|-------|
| `POST_MIN_GENES` | 500 | Raise if low-quality cells persist after filtering |
| `POST_MAX_MT` | 20.0 | Adjust after viewing QC violin plots |
| `POST_MAX_UMIS` | 100000 | Similarly look at QC violins to tune |
| `RUN_SCEVAN` | TRUE | Enable if **any** sample in the cohort is cancer/tumor — run on all samples including normals (they serve as the CNA reference baseline). Set FALSE only if the entire study is normal tissue. |
| `RUN_DECONTX` | TRUE | Recommended; set FALSE to skip ambient RNA correction |

### Scripts 03 & 04

| Parameter | Where | What to change |
|-----------|-------|----------------|
| `CONDITION_COLUMN` | Script 03, Section 1.3 | Metadata column for the two groups (e.g., `"Condition"`) |
| `GROUP_1` / `GROUP_2` | Script 03, Section 1.3 | `GROUP_1` = KO/treatment, `GROUP_2` = WT/control — positive log2FC = up in GROUP_1 |
| `CONDITION_LEVELS` | Script 04, Section 1.3 | Plot axis order: **WT/control first**, then KO/treatment — `c("WT", "KO")` |
| `TARGET_CELL_TYPE` | Script 04, Section 1.2 | Cell type to subset for scoring (e.g., `"T cells"`) |

> **Convention — DE direction vs. plot order:**
> `GROUP_1 = KO` in Script 03 so that positive log2FC always means *upregulated in KO*.
> `CONDITION_LEVELS = c("WT", "KO")` in Script 04 so that WT is always the visual baseline (left axis).
> These two settings are intentionally independent — changing one does not affect the other.

---

## Metadata File Format

The metadata Excel file is the single source of truth. Every column is
automatically attached to every cell from that sample.

**Required column:**

| Column | Requirement |
|--------|-------------|
| `SampleID` | Must exactly match the subfolder name inside `H5_DIR` |

**Recommended columns** (add any combination):

```
SampleID | Condition | Genotype | Diet | Sex | Timepoint | BatchID
S01_WT_M |    WT     |    WT    | ctrl |  M  |   Week4   |  Batch1
S02_KO_M |    KO     | Nr4a1_KO | ctrl |  M  |   Week4   |  Batch1
```

**Auto-created columns:**
- `group` = `Condition_Sex` (if both columns exist)
- Any `ColA_ColB` pattern → auto-built by Scripts 02–04 when referenced

---

## Marker Gene CSV Format

`cell_type_markers.csv` — one row per cell type/marker combination:

| Column | Description |
|--------|-------------|
| `cell_type` | Display name (e.g., `T cells`) |
| `parent_cell_type` | Broad type this sub-type belongs to (same as `cell_type` for broad markers) |
| `tier` | `broad` or `sub` |
| `markers` | Pipe-separated gene symbols: `Cd3e\|Cd3d\|Cd3g` |
| `sub_resolution` | Leiden resolution for sub-clustering (used in Script 02) |

---

## Mode D Enrichment — Key Knobs

All in `03_differential_expression.R` Part 1, Sections 1.8–1.9:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `GSEA_CATEGORY` / `GSEA_SUBCATEGORY` | `"C5"` / `"BP"` | MSigDB collection for GSEA (`"H"/""`= Hallmark, `"C2"/"KEGG"` = KEGG) |
| `GSEA_SPECIES` | `"Mus musculus"` | Change to `"Homo sapiens"` for human data |
| `GSEA_NPERM` | 1000 | Increase to 10,000 for publication figures |
| `ENRICHR_DATABASES` | GO BP 2025, GO MF 2025, KEGG 2026 | User-definable — verify names with `enrichR::listEnrichrDbs()` |
| `ENRICHR_TOP_N_GENES` | 250 | Top N significant DE genes (padj < 0.05) submitted to Enrichr per list; DV list uses no cutoff |

---

## Output Folder Structure

```
seurat_output/
├── <PROJECT>_processed_for_annotation.rds   # Script 01 → 02 handoff
├── <PROJECT>_final_annotated.rds            # Script 02 → 03/04 handoff
├── scevan_per_sample_results/               # Per-sample SCEVAN + checkpoints
├── QC_plots/                                # Pre/post-filter violin plots
├── annotation_plots/                        # UMAPs, dotplots, scoring heatmaps
├── differential_expression/
│   ├── ModeA_Wilcoxon_<G1>_vs_<G2>/        # DE results + volcano plots
│   └── ModeD_SplineDV_<G1>_vs_<G2>/
│       └── <CellType>/
│           ├── *_SplineDV_full_results.csv
│           ├── *_genes_DE_up_top250.csv
│           ├── *_genes_DE_dn_top250.csv
│           ├── *_genes_DV_top250.csv
│           ├── *_genes_overlap_DV_DEup.csv
│           ├── *_genes_overlap_DV_DEdn.csv
│           ├── *_GSEA_results.csv
│           ├── *_GSEA_NES_barplot.png
│           └── *_Enrichr_<DB>.png  (× 5 lists × 3 databases)
└── cell_scoring/
    ├── AUCell_scores_*.xlsx
    ├── AUCell_pseudobulk_ttest_*.xlsx
    └── pathway_plots/
```

---

## Common Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Script 01 crashes with OOM on > 10 samples | Memory exhausted during SCEVAN | Checkpoints allow resuming — restart after freeing RAM; SCEVAN runs one sample at a time |
| `DoubletFinder` removed > 14% of cells | Rollback triggered | Rollback is intentional — check pK sweep plot and increase `DF_PK_RANGE_MAX` if needed |
| Mode D Enrichr returns no results | Internet not available, or wrong DB name | Check connection; run `enrichR::listEnrichrDbs()` to verify database names |
| GSEA shows 0 significant pathways | Too few ranked genes or wrong species | Confirm `GSEA_SPECIES` matches your data; try `GSEA_SUBCATEGORY <- ""` for Hallmark |
| Annotation map gives wrong cell types | Weighted pre-scores misled annotation | Cross-check with DotPlot — pre-scores are a guide, not ground truth |

---

*Full parameter reference: Single-cell Standard Operating Procedure.docx — Sections 6–8*
