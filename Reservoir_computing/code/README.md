# scReservoir

[![PyPI version](https://img.shields.io/pypi/v/scReservoir)](https://pypi.org/project/scReservoir/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/screservoir/screservoir/blob/main/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1038/placeholder-blue)](https://doi.org/10.1038/placeholder)

**scReservoir** is a Python framework for inferring gene regulatory networks (GRNs) from single-cell RNA-seq data using reservoir computing (Echo State Networks). The method identifies attractor states, models developmental landscapes, and enables causal GRN inference through pseudotime-aware analysis.

## Key Features

- **Reservoir Computing**: Efficient nonlinear feature extraction using random recurrent neural networks
- **GRN Inference**: Three inference modes (standard, causal, velocity-based) for regulatory network discovery
- **Attractor Landscape**: Identification of cell fate attractors and Waddington energy landscapes
- **Scalability**: Sparse matrices and randomized SVD for computational efficiency on large datasets
- **Visualization**: Publication-ready plots for networks, landscapes, and gene programs
- **Pseudotime Integration**: Causal inference through temporal ordering of cells

## Installation

### From PyPI

```bash
pip install scReservoir
```

### From Source

```bash
git clone https://github.com/screservoir/screservoir.git
cd scReservoir
pip install -e .
```

## Quick Start

### Basic GRN Inference

```python
import scReservoir as scr
import scanpy as sc

# Load and preprocess single-cell data
adata = sc.datasets.pbmc68k_reduced()
X = scr.preprocess(adata, n_hvg=2000)

# Initialize reservoir computing model
reservoir = scr.ScReservoir(
    n_reservoir=500,
    spectral_radius=0.9,
    input_scale=0.5,
    leak_rate=0.3,
    density=0.01
)

# Compute reservoir states
H = reservoir.compute_states(X)

# Infer gene regulatory network
grn_model = scr.ScReservoirGRN(reservoir, lambda_reg=1e-3)
grn_model.fit(X, H, mode='standard')

# Get results
GRN = grn_model.get_grn()
top_regulators = grn_model.get_top_regulators(adata.var_names, k=10)

# Visualize
scr.plot_grn_heatmap(GRN, adata.var_names, top_n=50)
```

### Causal GRN with Pseudotime

```python
# Compute or load pseudotime
sc.tl.dpt(adata)
pseudotime = adata.obs['dpt_pseudotime'].values

# Sort data by pseudotime
X_sorted, pt_sorted = scr.order_by_pseudotime(X, pseudotime)

# Compute reservoir states
H = reservoir.compute_states(X_sorted, pseudotime=pt_sorted)

# Fit causal GRN (time-lagged)
grn_causal = scr.ScReservoirGRN(reservoir, lambda_reg=1e-3)
grn_causal.fit(X_sorted, H, pseudotime=pt_sorted, mode='causal')

GRN_causal = grn_causal.get_grn()
```

### Attractor Landscape and Fate Prediction

```python
# Initialize landscape model
landscape = scr.ScReservoirLandscape(
    n_attractors=5,
    svd_rank=50,
    lambda_reg=1e-3
)

# Fit landscape
landscape.fit(X_sorted, H, pseudotime=pt_sorted)

# Get results
energy = landscape.get_energy_landscape()
fate_probs = landscape.get_fate_probabilities()
attractors = landscape.get_attractor_states()
programs = landscape.get_attractor_gene_programs(X_sorted)

# Visualize
scr.plot_energy_landscape(energy, X_pca[:, :3], pseudotime=pt_sorted)
scr.plot_fate_probabilities(fate_probs, X_pca[:, :2])
scr.plot_attractor_genes(programs, adata.var_names, top_n=20)
```

## API Overview

### Core Classes

#### `ScReservoir`
Implements Echo State Network for high-dimensional nonlinear feature extraction.

**Key Methods:**
- `compute_states(X, pseudotime, washout)` - Run reservoir dynamics
- `compute_states_graph(X, adjacency, n_iter)` - Graph-propagated states

#### `ScReservoirGRN`
Gene regulatory network inference through regression on reservoir states.

**Key Methods:**
- `fit(X, H, pseudotime, mode)` - Fit GRN (modes: 'standard', 'causal', 'velocity')
- `get_grn()` - Normalized GRN matrix
- `get_top_regulators(gene_names, k)` - Top k regulators per gene
- `get_grn_sparse(threshold)` - Thresholded GRN

#### `ScReservoirLandscape`
Attractor detection and developmental landscape modeling.

**Key Methods:**
- `fit(X, H, pseudotime, velocity)` - Fit landscape model
- `get_energy_landscape()` - Energy values per cell
- `get_fate_probabilities()` - Fate commitment probabilities
- `get_attractor_states()` - Attractor gene expression profiles
- `get_grn_from_dynamics(W_in)` - GRN derived from latent dynamics

#### `ScGraphReservoir`
Graph-aware reservoir computing using k-NN cell neighborhoods.

### Utility Functions

- `preprocess(adata, n_hvg, log_transform)` - Data preprocessing
- `order_by_pseudotime(X, pseudotime)` - Sort cells by pseudotime
- `compute_gene_velocity(X, pseudotime, window)` - dX/dt estimation
- `ridge_regression(H, y, lambda_reg)` - Regularized regression
- `normalize_grn(GRN, method)` - GRN normalization
- `threshold_grn(GRN, threshold)` - Sparsification
- `compute_sensitivity(GRN)` - In-degree analysis
- `compute_activity(GRN)` - Out-degree analysis

### Visualization Functions

- `plot_grn_heatmap(GRN, gene_names, top_n, cmap)` - Network heatmap
- `plot_energy_landscape(energy, embedding, pseudotime)` - Energy landscape
- `plot_fate_probabilities(fate_probs, embedding, attractor_names)` - Fate probabilities
- `plot_attractor_genes(programs, gene_names, top_n)` - Attractor programs
- `plot_reservoir_dynamics(H, pseudotime, n_dims)` - Latent dynamics

## Method Description

scReservoir adapts **Echo State Networks** to single-cell transcriptomics:

1. **Random Reservoir**: A sparsely-connected recurrent neural network with fixed random weights. Only the linear readout is trained, reducing computational cost.

2. **GRN Inference**: For each gene, regression coefficients from reservoir states are back-projected through input connections to infer gene influence relationships: `GRN[i,j] ≈ |W_in[i] · W_out[:,j]|`

3. **Causal GRN**: Pseudotime-based ordering enables time-lagged regression, inferring directional regulatory relationships from temporal dynamics.

4. **Velocity Integration**: Models gene expression dynamics (dX/dt) instead of snapshots, improving causality inference.

5. **Attractor Landscape**:
   - Compressed latent dynamics via randomized SVD
   - Low-velocity cells identified as attractors
   - Energy landscape computed via kernel density estimation
   - Fate probabilities from distances to attractor centers

6. **Scalability**:
   - Sparse reservoir weights for memory efficiency
   - Randomized SVD for rank reduction
   - Compatible with large-scale datasets

## Examples

Three complete examples are provided in `examples/`:

1. **01_basic_grn_inference.py** - Standard GRN inference on synthetic data
2. **02_pseudotime_causal_grn.py** - Causal GRN with pseudotime ordering
3. **03_attractor_landscape.py** - Developmental landscape and fate prediction

Run examples:

```bash
cd examples
python 01_basic_grn_inference.py
python 02_pseudotime_causal_grn.py
python 03_attractor_landscape.py
```

## Citation

If you use scReservoir in your research, please cite:

```bibtex
@article{screservoir2024,
  title={scReservoir: Reservoir Computing for Single-Cell Gene Regulatory Network Inference},
  author={Contributors, scReservoir},
  journal={Nature Methods},
  year={2024},
  doi={10.1038/placeholder}
}
```

## Requirements

- Python >= 3.8
- numpy >= 1.21
- scipy >= 1.7
- scikit-learn >= 1.0
- scanpy >= 1.9
- anndata >= 0.8
- matplotlib >= 3.5
- seaborn >= 0.11

## License

scReservoir is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please open issues and pull requests on [GitHub](https://github.com/screservoir/screservoir).

## Support

For questions and discussions:
- **Issues**: [GitHub Issues](https://github.com/screservoir/screservoir/issues)
- **Documentation**: [ReadTheDocs](https://screservoir.readthedocs.io)
- **Email**: contact@screservoir.org

## References

- **Echo State Networks**: Jaeger et al., "The "echo state" approach to analysing and training recurrent neural networks" (2001)
- **Waddington Landscape**: Waddington, "The Strategy of the Genes" (1957)
- **Reservoir Computing**: Lukoševičius & Jaeger, "Reservoir Computing Approaches for Representation and Learning in Boolean Functions" (2009)

## Changelog

### Version 0.1.0 (Initial Release)
- Core reservoir computing implementation
- Standard, causal, and velocity-based GRN inference
- Attractor landscape modeling
- Graph-aware reservoir computing
- Comprehensive visualization tools
- Three complete examples

## Authors and Acknowledgments

scReservoir is developed by the Computational Biology and Systems Research community. We gratefully acknowledge contributions from [list contributors].

---

For the latest version and updates, visit: https://github.com/screservoir/screservoir
