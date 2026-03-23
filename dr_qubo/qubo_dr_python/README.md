# QUBO-Based Differential RNA Co-expression Analysis (`qubo_dr`)

A Python implementation of a quantum-inspired optimization pipeline for identifying gene subnetworks with maximum differential co-expression between biological conditions (e.g., Knockout vs. Wild Type) within specified biological pathways.

## Overview

This package casts the subnetwork selection problem as a **Quadratic Unconstrained Binary Optimization (QUBO)** problem, solvable via classical simulated annealing or quantum annealers (D-Wave). The algorithm is based on computational biology principles for scRNA-seq data analysis.

### Key Features

- **Pathway Gene Fetching**: Automatic retrieval from KEGG and Gene Ontology (GO Biological Process) databases
- **Differential Co-expression**: Gram matrix-based similarity computation across conditions
- **Graph-Based Adjacency**: Mutual Nearest Neighbor (MNN) or KNN sparse adjacency construction
- **QUBO Formulation**: Cardinality-constrained optimization with pathway network penalties
- **Flexible Solving**: Simulated annealing with graceful fallback (no D-Wave required)
- **Publication-Quality Visualization**: Co-expression heatmaps and force-directed network graphs

## Installation

### Basic Setup

```bash
pip install -r requirements.txt
```

### Optional Dependencies

- **D-Wave Quantum Tools** (for advanced annealing): `pip install dimod dwave-neal`
- **AnnData Support** (for .h5ad files): `pip install anndata`

If D-Wave tools are not installed, the pipeline will automatically use a pure NumPy fallback simulated annealing implementation.

## Quick Start

### Command Line Interface

```bash
# Run pipeline with KEGG Wnt signaling pathway (human)
python run_pipeline.py \
    --data counts.h5ad \
    --pathway-source kegg \
    --pathway-id 04310 \
    --organism hsa \
    --K 30 \
    --num-reads 2000

# Run with GO Biological Process
python run_pipeline.py \
    --data counts.h5ad \
    --pathway-source gobp \
    --pathway-id "Wnt signaling pathway" \
    --K 30

# List available KEGG pathways
python run_pipeline.py --list-kegg-pathways

# Search GO BP terms
python run_pipeline.py --list-gobp-terms wnt
```

### Python API

```python
from qubo_dr import run_pipeline
import anndata
import numpy as np

# Load scRNA-seq data
adata = anndata.read_h5ad('counts.h5ad')
X = adata.X.T.toarray()  # genes × cells
g = adata.var_names.tolist()
batch_id = adata.obs['batch'].values

# Run pipeline
results = run_pipeline(
    X=X,
    g=g,
    batch_id=batch_id,
    K=30,
    source='kegg',
    pathway_id='04310',  # Wnt signaling
    organism='hsa',
    ko_label='KO',
    wt_label='WT',
    plotit=True,
)

# Access results
selected_genes = results['sub_g_net']
subnetwork_matrix = results['sub_Q_net']
figures = results['figures']

print(f"Selected {len(selected_genes)} genes: {selected_genes}")
```

## Algorithm Overview

### Input
- Normalized scRNA-seq count matrix X (genes × cells)
- Gene names g and condition labels per cell
- Target pathway gene list (from KEGG/GO or user-provided)
- Target subnetwork size K

### Processing Steps

1. **Library-Size Normalization**: Per-cell scaling and log transformation
2. **Condition Subsetting**: Extract KO and WT cell submatrices
3. **Similarity Matrices**: Row-normalized co-expression (cosine similarity via Gram matrices)
4. **Differential**: Xdiff = Xwt_coexpr - Xko_coexpr (diagonal = 0)
5. **Graph Construction**: SVD → KNN → MNN adjacency for each condition
6. **Pathway Network**: Sparse matrix encoding pathway gene pairs (-1 penalty)
7. **QUBO Assembly**: Q1 = Xdiff + diag(Vdiff) - (MNN_wt + MNN_ko) + Xnet_target
   - Add cardinality penalty: Q = Q1 + P*(off-diagonal) + P*(1-2K)*diagonal
8. **Solving**: Minimize E(z) = z^T Q z via simulated annealing (z ∈ {0,1}^n)
9. **Extraction**: Subnetwork co-expression matrix and network graph
10. **Visualization**: Heatmaps and force-directed network with condition-specific edge colors

### Output
- Selected gene names (text file)
- Subnetwork co-expression matrix
- Network graph (networkx.Graph)
- Figures (heatmaps, network graphs)

## Package Structure

```
qubo_dr/
├── __init__.py          # Package initialization, version
├── preprocess.py        # Data normalization, similarity matrices
├── graph.py             # KNN/MNN adjacency construction
├── pathway.py           # KEGG & GO database interface
├── qubo.py              # QUBO matrix assembly and solving
├── plot.py              # Visualization functions
└── pipeline.py          # End-to-end orchestration
run_pipeline.py          # CLI entry point
requirements.txt         # Python dependencies
```

## Module Reference

### `qubo_dr.preprocess`
- `normalize_libsize()`: Library-size normalization + log transformation
- `subset_by_condition()`: Extract cells matching a condition
- `compute_gene_similarity()`: Gram matrix (cosine similarity)
- `compute_differential()`: Differential co-expression matrices

### `qubo_dr.graph`
- `build_mnn_adjacency()`: Construct MNN or KNN sparse adjacency matrix

### `qubo_dr.pathway`
- `fetch_kegg_pathway()`: Retrieve gene list from KEGG REST API
- `fetch_kegg_pathway_list()`: List all KEGG pathways for organism
- `search_kegg_pathways()`: Keyword search in KEGG pathways
- `fetch_gobp_geneset()`: Retrieve GO BP genes (gseapy or QuickGO API)
- `list_gobp_terms()`: Search GO BP terms by keyword
- `get_pathway_genes()`: Unified interface for KEGG/GOBP

### `qubo_dr.qubo`
- `build_pathway_network_matrix()`: Encode pathway structure as matrix
- `assemble_qubo_matrix()`: Construct Q with cardinality constraint
- `solve_qubo_simulated_annealing()`: Minimize QUBO via SA
- `extract_subnetwork()`: Select matrices for chosen genes

### `qubo_dr.plot`
- `plot_coexpr_heatmap()`: Co-expression heatmap with divergent colormap
- `plot_gene_network()`: Force-directed network with colored edges
  - Blue edges: WT co-expression gain
  - Red edges: KO co-expression gain
- `plot_condition_heatmaps()`: 3-panel comparison (KO, WT, Differential)

### `qubo_dr.pipeline`
- `run_pipeline()`: End-to-end analysis orchestration

## Common Use Cases

### 1. KEGG Pathway Analysis (Human)

```python
results = run_pipeline(
    X=X, g=g, batch_id=batch_id,
    source='kegg', pathway_id='04310',  # Wnt signaling
    organism='hsa',
    K=30
)
```

### 2. GO Biological Process Analysis

```python
results = run_pipeline(
    X=X, g=g, batch_id=batch_id,
    source='gobp', pathway_id='Wnt signaling pathway',
    organism='human',
    K=25
)
```

### 3. Custom Gene List

```python
my_genes = ['WNT1', 'WNT3', 'GSK3B', 'CTNNB1', 'APC', ...]

results = run_pipeline(
    X=X, g=g, batch_id=batch_id,
    genelist=my_genes,
    K=20
)
```

### 4. Cell State Integration

```python
# Include cell state trajectory information
results = run_pipeline(
    X=X, g=g, batch_id=batch_id,
    source='kegg', pathway_id='04310',
    use_cell_state=True,
    cell_state=pseudotime_values,  # Cells in trajectory order
    K=30
)
```

## Command-Line Arguments Reference

```
Pipeline Arguments:
  --data DATA                    Path to .h5ad or .npy file
  --pathway-source {kegg,gobp}   Database for pathway
  --pathway-id ID                KEGG ID (e.g., 04310) or GO term name
  --organism ORGANISM            Organism for mapping (default: human/hsa)
  --K K                          Subnetwork size (default: 30)
  --ko-label LABEL              KO condition in batch column (default: KO)
  --wt-label LABEL              WT condition in batch column (default: WT)
  --method {knn,mnn}            Graph construction (default: mnn)
  --neighbors N                  Number of neighbors (default: 30)
  --num-reads N                  SA iterations (default: 1000)
  --no-plot                     Suppress plot generation
  --outfile FILE                Output filename (default: qubo_genes_solution.txt)

Utility Arguments:
  --list-kegg-pathways          List KEGG pathways and exit
  --list-gobp-terms KEYWORD     Search GO BP terms and exit
```

## Input Data Format

### AnnData (.h5ad)
- `X`: count matrix (cells × genes, typically stored as sparse)
- `var_names`: gene names
- `obs['batch']`: condition labels per cell (e.g., 'KO', 'WT')

### NumPy (.npy)
- Store genes × cells matrix directly
- Use Python API with separate `g` and `batch_id` arrays

## Output Files

### Default: `qubo_genes_solution.txt`
```
# QUBO differential co-expression analysis results
# Selected genes: 30
#
GENE1
GENE2
GENE3
...
```

## Advanced Usage

### Tuning QUBO Parameters

```python
results = run_pipeline(
    ...,
    K=40,              # Larger subnetwork
    method='mnn',      # Use mutual nearest neighbors (stricter)
    n_neighbors=20,    # Fewer neighbors per gene
    n_svd=100,         # More SVD components
    penalty_scale=15,  # Stronger cardinality penalty
    num_reads=5000,    # More SA annealing iterations
)
```

### Network Analysis with NetworkX

```python
import networkx as nx
import matplotlib.pyplot as plt

G = results['G_graph']

# Node importance (from QUBO objective)
node_sizes = 300 * (results['sub_Qv'] / results['sub_Qv'].max())

# Compute centrality
betweenness = nx.betweenness_centrality(G)
print("Top genes by betweenness centrality:")
for gene, score in sorted(betweenness.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {results['sub_g_net'][gene]}: {score:.3f}")
```

### Visualization Customization

```python
from qubo_dr.plot import plot_coexpr_heatmap, plot_gene_network
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 6))

# Custom heatmap
plot_coexpr_heatmap(
    results['sub_Q_net'],
    results['sub_g_net'],
    title='My Custom Title',
    cmap='coolwarm',
    ax=ax
)

fig.savefig('heatmap.png', dpi=300, bbox_inches='tight')
```

## Troubleshooting

### "No genes found matching pathway"
- Verify pathway ID format (KEGG: bare ID like '04310' or prefixed 'hsa04310')
- Check organism code (KEGG: 'hsa', 'mmu', 'dme'; GOBP: 'human', 'mouse')
- Ensure gene names in data match database format (often case-sensitive)

### "No cells found matching condition label"
- Verify condition labels in batch column match `--ko-label` and `--wt-label`
- Check unique values: `adata.obs['batch'].unique()`

### ImportError: dimod
- Optional dependency. Pipeline uses NumPy fallback automatically.
- For D-Wave support: `pip install dimod dwave-neal`

### ImportError: anndata
- Optional dependency. Use `--data file.npy` with manual loading, or `pip install anndata`

## Citation

If you use this package in research, please cite:

```
QUBO-based Differential Co-expression Analysis for scRNA-seq
Author: Selim Romero, Texas A&M University
```

## License

[Specify license here]

## Contact

Selim Romero
Texas A&M University
[Email/GitHub]

## Acknowledgments

- KEGG REST API: https://www.kegg.jp/kegg/rest/
- Gene Ontology & QuickGO: https://www.ebi.ac.uk/QuickGO/
- gseapy library: https://gseapy.readthedocs.io/
- D-Wave quantum annealing: https://www.dwavesys.com/

## References

- Differential co-expression analysis concepts
- QUBO formulation for optimization
- scRNA-seq preprocessing methodologies
- Graph-based gene network analysis
