# TimeSCape v0.2 — Analysis Summary

## What was built and why

---

### 1. Circadian gene identification (Section 5)

Cosinor regression was fitted to each gene's expression trajectory across ZT time points within the target cell type. Genes with a significant rhythmic fit (adjusted p-value below threshold) and a minimum amplitude were retained as **confident circadian genes**. This step answers the question: *which genes genuinely oscillate in these cells, and when do they peak?*

---

### 2. Phase-bin ORA — assigning pathways to time windows (Section 7–8)

Confident circadian genes were grouped into phase bins (default 1-hour windows) according to their cosinor-estimated acrophase. Within each bin, a **hypergeometric over-representation analysis (ORA)** tested whether any gene set — drawn from KEGG, Reactome, and GO Biological Process — was significantly enriched among genes peaking in that window.

An **EnrichR consensus layer** was optionally added: the same bin genes were queried against the EnrichR API (KEGG 2026, Reactome 2024, GO:BP 2025), and only pathways significant in **both** phyper and EnrichR were labelled "consensus", giving higher confidence. Pathways significant only in phyper were retained as "phyper only" — they are still informative but less certain.

**Why three databases?** KEGG provides compact, curated metabolic and signalling pathways. Reactome provides mechanistic detail on immune signalling, relevant for CD8⁺ T cells. GO:BP captures biological process terms that neither KEGG nor Reactome cover (e.g. RNA splicing regulation, protein folding rhythms). Using all three maximises coverage while the dual-ORA filter limits false positives.

**Why no exclusion patterns?** The downstream cosinor scoring (Section 10) naturally rejects pathways that do not show oscillatory AUC scores across ZT, so there is no need to pre-filter by name.

---

### 3. AUCell phase-restricted scoring (Section 9)

For each pathway that emerged from the ORA, only the subset of its genes that peak in the same phase bin were used to score cells (**phase-restricted gene sets**). This is the key innovation: conventional pathway scoring pools all members regardless of their peak time, which causes phase cancellation — a gene peaking at ZT6 and one peaking at ZT18 cancel each other's signal. By restricting to co-phased genes, the AUCell score captures the actual temporal activity of the pathway rather than a time-averaged blur.

---

### 4. Pathway cosinor — identifying oscillating pathways (Section 10)

AUCell scores per cell were treated as a new expression matrix and fitted with cosinor regression, exactly as was done for individual genes. This identified which pathways oscillate at the cell-population level and estimated their acrophase (time of peak activity) and amplitude. Violin plots across ZT time points with overlaid cosine fits allowed visual confirmation.

---

### 5. Gene Regulatory Network time series (Section 11)

For each confident oscillating pathway, a **gene regulatory network (GRN)** was constructed at every ZT time point using Pearson correlations between co-expressed genes. The network nodes are the union of:

- **Top 1–3 circadian representatives per phase bin** — genes with the strongest amplitude + mesor within each ZT window, ensuring every part of the circadian cycle is represented
- **Detected core clock genes** (Arntl/Bmal1, Clock, Per1/2/3, Cry1/2, Nr1d1/2, Rora/b/c, Dbp, Tef, Hlf, etc.) present among the confident circadian genes
- **Phase-restricted pathway genes** — the members of this specific pathway that peak in the relevant bin

By examining how edges appear and disappear across ZT00 → ZT21, you can see when the regulatory wiring of each pathway is most active, which hub genes gain or lose connections across the day, and whether clock genes are central drivers or peripheral modulators in each pathway context.

**Node colours:** orange-red = circadian oscillator; blue = pathway gene; purple = both; grey = co-expressed bystander.  
**Edge colours:** red = positive co-expression; blue = negative (anti-correlated).

---

### Output files

| Section | Output |
|---|---|
| 5 | Cosinor results table (`T1`), confident gene list |
| 7–8 | `phase_results` — ORA tables per bin, consensus labels |
| 9 | AUCell matrix (cached as `.rds`), phase-restricted gene sets |
| 10 | Pathway cosinor stats, violin plots per pathway |
| 11 | One GRN time-series PNG per confident pathway, saved to `<cell_type>_GRN/` |

---

*TimeSCape v0.2 — developed for the Romero Lab, Texas A&M University*
