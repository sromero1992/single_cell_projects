# QUBO Differential Co-expression — Pathway Subnetwork Discovery in scRNA-seq

**Author:** Selim Romero, Texas A&M University

---

## Overview

This pipeline identifies the subset of genes within any biological pathway that carries the **maximum differential co-expression signal** between two cellular conditions. Given two populations of cells — such as knockout vs. wild type, disease vs. healthy, treated vs. control, or two distinct cell states — the method finds the core gene subnetwork whose co-regulatory relationships have shifted most between conditions.

The selection problem is cast as a **Quadratic Unconstrained Binary Optimization (QUBO)**, making it directly compatible with D-Wave quantum annealers and solvable classically via simulated annealing.

---

## The Core Idea

Standard differential expression asks: *which genes changed in average level?*

This pipeline asks something fundamentally different: *which genes changed in how they relate to each other?*

Two genes can maintain the same mean expression level across conditions while completely rewiring their co-regulatory relationship. Those rewired edges often represent the mechanistic core of a perturbation — the signaling bottleneck, the feedback loop that breaks, the hub that loses or gains influence. By searching for the subnetwork whose pairwise co-expression has shifted most, the method recovers biologically actionable structure that univariate differential expression misses.

The QUBO objective jointly optimizes four terms, plus a hard cardinality constraint:

1. **Pairwise differential co-expression** (`X_diff = S_B − S_A`) — which gene pairs gained co-expression in condition A (test) relative to condition B (reference). Computed as the difference of non-negative cosine similarity (Gram) matrices, so all values are bounded and directionality is preserved.

2. **Local co-expression topology (MNN)** — Mutual Nearest Neighbor adjacency matrices built from a truncated SVD embedding of each condition's expression matrix. The SVD is run on cells, giving each gene a coordinate in cell-space; KNN is then between genes. Subtracting these matrices rewards selecting genes that are local co-expression hubs within each condition — structurally connected, not just correlated on average.

3. **Cell-state difference vector `V_diff`** (optional linear bias) — for each gene, the difference in how strongly its expression aligns with a per-cell biological state scalar between the two conditions. The scalar can be **any continuous per-cell measure**: UCell pathway activity score, pseudotime coordinate, cell potency score, differentiation rank, or similar. Pass `None` / `[]` to omit.

4. **Subnetwork size constraint** — exact selection of exactly K genes is enforced via a quadratic cardinality penalty embedded directly in the QUBO matrix, requiring no external solver constraint handling.

**A note on the pathway membership prior (`X_net`):**
This optional term assigns a −1 energy bonus to co-selected gene pairs that are both core pathway members. When the expression matrix is already restricted to pathway genes (the standard case), every off-diagonal pair qualifies and the matrix becomes a uniform constant — it shifts the energy landscape equally for all solutions and has no effect on which solution is optimal. `X_net` is therefore **off by default** (`use_pathway_prior=False`). Enable it only when you append extra candidate genes (TF targets, GWAS hits, drug targets) to the pathway gene list — in that mixed case, it creates a preference for core-pathway pairs over candidate-gene pairs.

---

## Reference Dataset: Gut Cell Atlas (Pan-GI Extended+)

The primary dataset used to develop and demonstrate this pipeline is the **[Gut Cell Atlas pan-gastrointestinal (pan-GI) Extended+ resource](https://www.gutcellatlas.org/pangi.html)** (Oliver et al., *Nature* 2024).

| Property | Value |
|---|---|
| Total cells | ~1.6 million |
| Healthy reference | ~1.1 million cells, 189 donors, 136 cell states |
| Disease datasets | 12 studies — coeliac disease, ulcerative colitis, Crohn's disease, GI cancers |
| Tissues covered | Oral mucosa, oesophagus, stomach, small intestine, large intestine, mesenteric lymph nodes |
| Integrated studies | 25 scRNA-seq datasets uniformly processed with scAutoQC |
| Access | [gutcellatlas.org/pangi.html](https://www.gutcellatlas.org/pangi.html) |

**Why this dataset?** The pan-GI atlas spans the full gastrointestinal axis across both healthy and diseased states, making it ideal for testing pathway co-expression rewiring in context. For example, Wnt/β-catenin signaling in colon epithelial stem cells can be interrogated across healthy donors and IBD patients within the same unified cell-type framework.

**Citation:**
> Oliver AJ, Huang N, Bartolome-Casado R, et al. Single-cell integration reveals metaplasia in inflammatory gut diseases. *Nature*. 2024;635(8039):699–707. doi:[10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1). PMID: 39567783.

---

## Biological Contexts

The pipeline is general. You supply a pathway and a pair of conditions — the rest adapts automatically. Below are representative use cases, but any scRNA-seq or bulk RNA-seq dataset with two annotated conditions and a relevant pathway can be analyzed.

| Cell type / Tissue | Pathway of interest | Condition A (test) → B (reference) | KEGG / GO term |
|---|---|---|---|
| Colon epithelial stem cells | Wnt / β-catenin signaling | Nuclear receptor KO → WT | `KEGG:hsa04310` |
| Colon stem cells (validation) | Wnt / β-catenin signaling | APC mutation → WT | `KEGG:mmu04310` |
| Cancer-associated fibroblasts (CAFs) | IGF / PI3K–Akt signaling | CAF vs. normal fibroblast | `KEGG:hsa04151` |
| Fibrotic lung fibroblasts | TGF-β / ECM remodeling | IPF vs. healthy | `KEGG:hsa04350` |
| Hepatic stellate cells | TGF-β / Hippo signaling | Activated vs. quiescent | `GO:0007179` |
| Macrophages / microglia | NF-κB / inflammatory | LPS-stimulated vs. resting | `KEGG:hsa04064` |
| T cells | TCR signaling / IL-2 | Exhausted vs. effector | `GO:0050852` |
| Pancreatic β-cells | Insulin secretion | Diabetic vs. healthy | `KEGG:hsa04911` |
| Neurons | Neurotrophin / MAPK | Alzheimer's vs. control | `KEGG:hsa04722` |

The pathway gene list functions as a **biologically informed prior**, not a hard filter. You can always augment any database pathway with additional candidate genes — transcription factors, drug targets, GWAS hits, or any gene you hypothesize plays a role — by simply appending them to the gene list before running.

---

## Pathway Sources

Pathways can come from three places, and all three can be combined.

### Option A — KEGG REST API (Python)

```python
from qubo_dr.pathway import get_pathway_genes, search_kegg_pathways

# Wnt signaling, mouse
genelist = get_pathway_genes('kegg', '04310', organism='mmu')

# PI3K-Akt (IGF / growth factor signaling), human
genelist = get_pathway_genes('kegg', '04151', organism='hsa')

# TGF-beta signaling
genelist = get_pathway_genes('kegg', '04350', organism='hsa')

# Explore available pathways
search_kegg_pathways('fibrosis', organism='hsa')
search_kegg_pathways('growth factor', organism='hsa')
search_kegg_pathways('Wnt', organism='mmu')
```

KEGG supports most model organisms. Use `'hsa'` for human, `'mmu'` for mouse, `'rno'` for rat, `'dme'` for fly, and so on.

### Option B — GO Biological Process (Python)

```python
from qubo_dr.pathway import get_pathway_genes, list_gobp_terms

# By term name (fuzzy matched)
genelist = get_pathway_genes('gobp', 'Wnt signaling pathway')
genelist = get_pathway_genes('gobp', 'fibroblast growth factor receptor signaling')
genelist = get_pathway_genes('gobp', 'TGF-beta receptor signaling pathway')
genelist = get_pathway_genes('gobp', 'insulin-like growth factor receptor signaling')

# Explore GO BP terms by keyword
list_gobp_terms('growth factor')
list_gobp_terms('fibrosis')
list_gobp_terms('stem cell')
```

### Option C — Custom list (MATLAB or Python)

```matlab
% MATLAB — define manually or load from a file
genelist = ["WNT5A" "FZD1" "CTNNB1" "DVL2" "GSK3B" "AXIN2" "LRP6" "TCF7L2" ...
            "LEF1"  "SFRP1" "TCF7" "FZD4" "FZD9" "BCL9" "LRP5"];
```

```python
# Python
genelist = ["WNT5A", "FZD1", "CTNNB1", "DVL2", "GSK3B", "AXIN2", "LRP6", "TCF7L2"]
```

### Adding Extra Genes to Any Pathway

You can always extend a database-derived pathway with additional candidates. This is useful when you have prior biological knowledge — a gene from a GWAS, a drug target, a co-factor known to interact with the pathway, or a transcription factor your lab is studying.

```matlab
% MATLAB — start from KEGG list, add candidates
genelist = ["WNT5A" "CTNNB1" "DVL2" ...];   % your base pathway
extra    = ["YOUR_TF" "YOUR_RECEPTOR" "CANDIDATE_GENE"];
genelist = [genelist, extra];
```

```python
# Python — same idea
base_genelist = get_pathway_genes('kegg', '04310', organism='hsa')
extra_genes   = ['YOUR_TF', 'YOUR_RECEPTOR', 'CANDIDATE_GENE']
genelist      = base_genelist + extra_genes
```

The optimizer will evaluate these extra genes in the context of the pathway — if they carry strong differential co-expression signal with pathway members, they will be selected; if not, they will be excluded.

---

## Input Data Format

The pipeline takes:

| Input | Type | Description |
|---|---|---|
| `X` | `G × N` numeric matrix | Normalized count matrix (genes × cells). Library-size normalize to ~10,000 counts per cell, optionally log1p-transform. |
| `g` | `1 × G` string array / list | Gene name vector. Names must be consistent with pathway gene names (case-insensitive matching is applied). |
| `batch_id` | `1 × N` string array / list | Condition label per cell. Any two-group labeling works: `'KO'`/`'WT'`, `'IPF'`/`'Healthy'`, `'Treated'`/`'Control'`, `'CAF'`/`'NF'`, etc. Partial matching is supported (e.g., `'KO'` will match `'Nr4a1_KO'`). |
| `genelist` | string array / list | Pathway gene list. Fetched from KEGG/GOBP or provided manually. Extra candidate genes can be appended. |
| `K` | integer | Desired subnetwork size — how many genes to select. Typically 15–50. |

---

## Usage

### MATLAB

```matlab
% Add package to path
addpath('/path/to/dr_qubo');

% Normalize your count matrix if not already done
X = pkg.norm_libsize(full(sce.X), 1e4);

% Fetch or define pathway gene list
genelist = ["WNT5A" "FZD1" "CTNNB1" "DVL2" "GSK3B" "AXIN2" "LRP6" "TCF7L2" ...
            "LEF1"  "SFRP1" "TCF7"   "FZD4" "FZD9"  "BCL9"  "LRP5" "FZD8"];
% Add any extra candidates
genelist = [genelist, "GENE_A", "GENE_B"];

% Configure options
opts.condition_A  = 'KO';       % label for condition A in batch_id
opts.condition_B  = 'WT';       % label for condition B
opts.K            = 30;         % subnetwork size
opts.n_neighbors  = 30;         % neighbors for MNN graph
opts.method       = 'mnn';      % 'mnn' (recommended) or 'knn'
opts.use_cell_state = false;    % set true to include cell state bias
opts.plotit       = true;       % generate heatmap and network figures
opts.outfile      = 'selected_genes.txt';

% Run
results = qubo_dr.run_pipeline(X, g, batch_id, genelist, opts);

% Inspect results
disp(results.sub_g_net)         % selected gene names
imagesc(results.sub_Q_net)      % differential co-expression submatrix
```

### Python

```python
from qubo_dr.pipeline import run_pipeline
from qubo_dr.pathway import get_pathway_genes

# ── Example 1: Wnt signaling, colon stem cells ─────────────────────────────
genelist = get_pathway_genes('kegg', '04310', organism='mmu')
# Optionally add genes of interest beyond the canonical pathway
genelist += ['Nr4a1', 'Ncoa3', 'Rarg']

results = run_pipeline(
    X          = X,           # (G, N) numpy array, normalized
    g          = g,           # list of G gene names
    batch_id   = batch_id,    # list of N condition labels
    genelist   = genelist,
    K          = 30,
    ko_label   = 'KO',
    wt_label   = 'WT',
    method     = 'mnn',
    n_neighbors = 30,
    plotit     = True,
    outfile    = 'selected_genes.txt'
)
print("Selected subnetwork:", results['sub_g_net'])


# ── Example 2: IGF / PI3K–Akt, cancer-associated fibroblasts ──────────────
genelist_igf = get_pathway_genes('kegg', '04151', organism='hsa')
genelist_igf += ['TGFB1', 'ACTA2', 'FAP']   # add fibroblast markers

results_caf = run_pipeline(
    X          = X_caf,
    g          = g_caf,
    batch_id   = condition_caf,
    genelist   = genelist_igf,
    K          = 25,
    ko_label   = 'CAF',
    wt_label   = 'NF',
    plotit     = True
)


# ── Example 3: TGF-β signaling, fibrotic lung fibroblasts ─────────────────
genelist_tgf = get_pathway_genes('gobp', 'TGF-beta receptor signaling pathway')

results_fib = run_pipeline(
    X          = X_fib,
    g          = g_fib,
    batch_id   = condition_fib,
    genelist   = genelist_tgf,
    K          = 20,
    ko_label   = 'IPF',
    wt_label   = 'Healthy',
    plotit     = True
)


# ── Example 4: NF-κB signaling, macrophage activation ─────────────────────
genelist_nfkb = get_pathway_genes('kegg', '04064', organism='hsa')

results_mac = run_pipeline(
    X          = X_mac,
    g          = g_mac,
    batch_id   = condition_mac,
    genelist   = genelist_nfkb,
    K          = 20,
    ko_label   = 'Stimulated',
    wt_label   = 'Resting',
    plotit     = True
)
```

### Command Line (Python)

```bash
# Install dependencies
cd qubo_dr_python && pip install -r requirements.txt

# Wnt pathway, mouse, from KEGG
python run_pipeline.py \
  --data mydata.h5ad \
  --pathway-source kegg \
  --pathway-id 04310 \
  --organism mmu \
  --K 30 \
  --condition-a KO --condition-b WT \
  --outfile selected_genes.txt

# IGF/PI3K pathway, human, from KEGG
python run_pipeline.py \
  --data fibroblast_data.h5ad \
  --pathway-source kegg \
  --pathway-id 04151 \
  --organism hsa \
  --K 25 \
  --condition-a CAF --condition-b NF

# TGF-beta signaling from GO Biological Process
python run_pipeline.py \
  --data lung_data.h5ad \
  --pathway-source gobp \
  --pathway-id "TGF-beta receptor signaling pathway" \
  --K 20 \
  --condition-a IPF --condition-b Healthy

# Explore what KEGG pathways are available (human)
python run_pipeline.py --list-kegg-pathways --organism hsa

# Search GO BP terms by keyword
python run_pipeline.py --list-gobp-terms "growth factor"
python run_pipeline.py --list-gobp-terms "fibroblast"
```

---

## Outputs

From the solution binary vector **z***, the pipeline produces:

| Output | Description |
|---|---|
| `sub_g_net` | The K selected gene names forming the most differentially co-expressed subnetwork |
| `sub_Q_net` | The differential co-expression matrix restricted to selected genes (G_selected × G_selected) |
| `sub_Qv` | Per-gene contribution scores — higher magnitude = more central to the subnetwork |
| Heatmap figure | Color-coded gene × gene differential co-expression for the subnetwork |
| Network figure | Force-directed graph; **red edges** = increased co-expression in condition A; **blue edges** = increased co-expression in condition B |
| `outfile` | Tab-separated text file of selected gene names |

The sign convention is: **negative values in `X_diff`** correspond to **increased co-expression in condition A** (first condition / KO). The optimizer finds genes whose pairwise co-expression is most negative — i.e., the genes most coordinately upregulated in condition A relative to condition B.

---

## Key Parameters

| Parameter | Description | Typical range | Notes |
|---|---|---|---|
| `K` | Subnetwork size | 15–50 | Start with ~30; reduce if the pathway is small |
| `n_neighbors` | Neighbors for MNN graph | 15–50 | For pathways < 100 genes, keep ≤ 30–40 to avoid a near-complete graph |
| `method` | Adjacency type | `'mnn'` / `'knn'` | MNN is stricter (mutual) and recommended |
| `penalty_scale` | Cardinality constraint weight | 5–20 | Multiplier on max\|Q1\|; higher = stricter K enforcement |
| `n_svd` | SVD components for dim. reduction | 30–100 | Reduce for very small gene sets |
| `use_cell_state` | Include cell state linear bias | `true` / `false` | Useful when comparing differentiation states |

---

## Statistical Validation

Running the QUBO gives you a subnetwork — but two questions remain: *Is this subnetwork significantly different from chance?* and *Are these specific genes reproducible?* Two procedures answer these:

### UCell Pathway Activity Scoring

Before and alongside the QUBO, each cell receives a pathway-specific activity score using the **UCell** method (Andreatta & Carmona, *Comput Struct Biotechnol J*, 2021). UCell ranks all genes by expression within each cell and computes a normalized Mann-Whitney U statistic for the pathway gene set. The result is a per-cell score in [0, 1] that is robust to zero-inflation and dataset size — no background gene set required.

This score serves two roles:
1. **Preliminary check**: a Mann-Whitney test on the UCell score distributions (condition A vs B) tells you whether the pathway is differentially active at all before investing in the QUBO. The AUROC between the two distributions gives an effect size.
2. **Linear QUBO bias (V_diff)**: the per-condition UCell score vectors are normalized and projected back onto gene space, yielding V_diff — a gene-level differential alignment vector that weights QUBO selection toward genes whose expression co-varies most with overall pathway activity.

This is preferred over general potency scores (like sc_potency) because it is pathway-specific: the linear bias reflects how much each gene drives the pathway's activity difference, rather than a cell's global transcriptional state.

```python
from qubo_dr.cell_scoring import compute_differential_activity

scores_A, scores_B, Vdiff, pval, auroc = compute_differential_activity(
    X, g, batch_id, genelist,
    condition_A='KO', condition_B='WT'
)
print(f"Pathway activity p-value (KO vs WT): {pval:.4f}")
print(f"AUROC: {auroc:.3f}")  # > 0.5 = higher activity in WT
```

```matlab
% MATLAB
scores = qubo_dr.stats.compute_ucell_score(X, g, genelist);
% scores: 1×N per-cell pathway activity in [0,1]
```

### Permutation Test (Subnetwork Significance)

The **label permutation test** answers: *is the selected subnetwork's differential co-expression stronger than expected under random condition assignment?*

The real QUBO solution z* is held fixed. Condition labels are shuffled randomly N_perm times; each time the differential co-expression matrix is recomputed and the real z* is scored against it. This builds a null distribution of subnetwork energies under the exchangeability null, requiring no QUBO re-solve per permutation.

- **p-value**: fraction of null energies as extreme as the real
- **z-score**: how many standard deviations below the null mean the real energy sits
- Rule of thumb: z < −3, p < 0.05 = the subnetwork captures a genuine condition-specific signal

```python
from qubo_dr.permutation import permutation_test

perm_results = permutation_test(
    X, g, batch_id,
    selected_idx = results['selected_idx'],
    Q1           = results['Q1'],
    n_perm       = 1000
)
print(f"p = {perm_results['pval']:.4f}, z = {perm_results['zscore']:.2f}")
```

```matlab
% MATLAB
opts.n_perm = 1000;
perm = qubo_dr.stats.permutation_test(X, g, batch_id, selected_idx, Q1, opts);
fprintf('p = %.4f,  z = %.3f\n', perm.pval, perm.zscore);
```

### Bootstrap Stability (Gene Selection Confidence)

The **bootstrap** answers: *which specific genes are robustly selected regardless of which cells happen to be in the sample?*

Cells are resampled with replacement independently within each condition (B = 200 iterations by default). The full QUBO pipeline is re-run on each bootstrap resample. The result is a per-gene **selection frequency** f_i ∈ [0, 1] — how often each gene is selected across bootstrap datasets.

- **Stable genes** (f_i > 0.5): selected in the majority of bootstrap samples → high confidence
- **Stability index Φ**: mean frequency of stable genes → overall subnetwork reproducibility (Φ near 1.0 = very stable)
- Genes in the real solution but with low bootstrap frequency may represent structurally ambiguous or outlier-driven selections

```python
from qubo_dr.permutation import bootstrap_stability

boot_results = bootstrap_stability(
    X, g, batch_id, genelist, K=30,
    n_boot=200, freq_thresh=0.5
)
print("Stable genes:", boot_results['stable_genes'])
print(f"Stability index Φ = {boot_results['stability_index']:.3f}")
```

```matlab
% MATLAB
opts.n_boot      = 200;
opts.freq_thresh = 0.5;
boot = qubo_dr.stats.bootstrap_stability(X, g, batch_id, genelist, K, opts);
fprintf('Stability index = %.3f | Stable genes: %d\n', ...
    boot.stability_index, numel(boot.stable_genes));
```

### Interpreting the Full Results Together

| Result | Interpretation |
|---|---|
| UCell p < 0.05 | The pathway is genuinely differentially active between conditions — QUBO is warranted |
| Permutation z < −3, p < 0.05 | The selected subnetwork captures a real condition-specific co-expression signal |
| Stability Φ > 0.8 | The specific gene selection is highly reproducible; confidence in biological interpretation |
| Stable genes ∩ real solution | The "core" high-confidence subnetwork — lead candidates for follow-up |
| Real solution \ stable genes | Secondary genes — may reflect condition-specific outliers; interpret with caution |

---

## Methodology Summary

A full mathematical derivation is provided in `methodology.tex` (LaTeX, Overleaf-ready). The abbreviated pipeline is:

1. **Normalize** the count matrix to library size (CPM × 10⁴, optional log1p)
2. **Subset** to pathway genes; split cells into two condition-specific matrices X_A and X_B
3. **Cosine similarity**: Row-normalize each gene's expression vector (L₂), compute Gram matrix S = Z·Zᵀ ∈ [0,1]^(G×G) per condition
4. **Differential matrix**: X_diff = S_B − S_A (diagonal = 0; negative entry = increased co-expression in condition A)
5. **MNN adjacency**: SVD reduction → KNN search → mutual neighbor filtering, independently per condition
6. **QUBO assembly**: Q₁ = X_diff + diag(V_diff) − (MNN_A + MNN_B) + X_net\_target
7. **Penalty**: P = 10·max|Q₁|; Q(i≠j) += P; Q(i,i) += P·(1−2K)
8. **Solve**: Minimize E(z) = zᵀQz over z ∈ {0,1}^G
9. **Extract**: Selected genes, subnetwork matrix, visualization

---

## Repository Structure

```
dr_qubo/
├── +qubo_dr/                          ← MATLAB package
│   ├── run_pipeline.m                 ← main entry point
│   ├── +preprocess/
│   │   ├── subset_by_condition.m      ← split X by condition label
│   │   ├── compute_gene_similarity.m  ← L2 normalization + Gram matrix
│   │   └── compute_differential.m     ← X_diff = S_B − S_A, optional V_diff
│   ├── +graph/
│   │   ├── build_adjacency_sparse.m   ← SVD + KNN/MNN for full gene set
│   │   ├── build_adjacency_subset.m   ← SVD + KNN/MNN for pathway subset
│   │   └── build_pathway_network.m    ← X_net_target prior matrix
│   ├── +qubo/
│   │   ├── assemble_qubo_matrix.m     ← build Q1, apply cardinality penalty
│   │   ├── solve_qubo_problem.m       ← MATLAB qubo() + solve()
│   │   └── extract_subnetwork.m       ← extract genes, sub_Q_net, sub_Qv
│   └── +plot/
│       ├── plot_coexpr_heatmap.m      ← heatmap of selected submatrix
│       └── plot_gene_network.m        ← force-layout network, colored edges
│
├── qubo_dr_python/                    ← Python package
│   ├── qubo_dr/
│   │   ├── __init__.py
│   │   ├── preprocess.py              ← normalize, subset, similarity, differential
│   │   ├── graph.py                   ← SVD + KNN/MNN adjacency (scipy sparse)
│   │   ├── pathway.py                 ← KEGG REST API + GO-BP via gseapy/QuickGO
│   │   ├── qubo.py                    ← Q1 assembly, penalty, SA solver (dimod/neal)
│   │   ├── plot.py                    ← heatmap, network, 3-panel condition comparison
│   │   └── pipeline.py                ← run_pipeline() end-to-end function
│   ├── run_pipeline.py                ← CLI entry point (argparse)
│   └── requirements.txt
│
├── README.md                          ← this file
├── methodology.tex                    ← LaTeX source (Overleaf-ready)
└── qubo_genes_sol_dr.txt              ← example output (Wnt test case, stem cells)
```

---

## Test Case

The included `qubo_genes_sol_dr.txt` was generated using the **canonical Wnt signaling pathway** (`KEGG mmu04310`) applied to mouse colon stem cell scRNA-seq data, comparing a nuclear receptor knockout to wild type. It serves as a reference output for verifying the pipeline is functioning correctly on your system.

The selected 32-gene subnetwork recovers Frizzled receptors, Wnt ligands, scaffold proteins, secreted modulators, transcription factors, and co-activators — a biologically coherent result consistent with the known co-regulatory structure of the Wnt pathway. Colon stem cell Wnt signaling is well-characterized in the literature, making it a reliable benchmark context.

---

## Dependencies

### MATLAB
- MATLAB R2023a or later
- Statistics and Machine Learning Toolbox (`knnsearch`, `normr`)
- Optimization Toolbox (`qubo`, `solve`)
- A single-cell toolkit providing `norm_libsize` and optionally `sc_potency` (e.g., [scGEAToolbox](https://github.com/jamesjcai/scGEAToolbox) or a lab `+pkg` package)

### Python
```
numpy >= 1.24          scipy >= 1.10         scikit-learn >= 1.2
matplotlib >= 3.7      seaborn >= 0.12       networkx >= 3.0
requests >= 2.28       gseapy >= 1.0         pandas >= 2.0
dimod >= 0.12          dwave-neal >= 0.6     anndata >= 0.9
tqdm >= 4.64                                               # progress bars (bootstrap)
```

Install with: `pip install -r qubo_dr_python/requirements.txt`

---

## Citation

Developed by **Selim Romero**, Texas A&M University.
Citation information will be updated upon publication.
