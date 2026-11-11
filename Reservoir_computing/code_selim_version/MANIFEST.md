# scReservoir — Repository Manifest
### Version 0.2.0 · Classical + Quantum Reservoir Computing for scRNA-seq

---

## Directory Structure

```
code/
│
├── classical/                         Classical Echo State Network pipeline
│   ├── __init__.py          (34 ln)   Docstring — usage guide & workflow
│   ├── utils.py            (324 ln)   Data preprocessing & linear algebra helpers
│   ├── reservoir.py        (277 ln)   Build reservoir weights, compute H matrix
│   ├── grn.py              (252 ln)   GRN inference (3 regression modes)
│   ├── landscape.py        (403 ln)   Attractor landscape, energy, fate probs
│   └── plotting.py         (321 ln)   All visualization functions
│
├── quantum/                           Quantum Reservoir Computing pipeline
│   ├── __init__.py          (41 ln)   Docstring — both usage patterns explained
│   └── quantum_reservoir.py(870 ln)   Full quantum RC + NG-RC implementation
│
├── examples/                          Runnable end-to-end scripts
│   ├── 01_classical_grn.py (111 ln)   Classical ESN → GRN inference
│   ├── 02_classical_landscape.py(99 ln) Classical ESN → attractor landscape
│   └── 03_quantum_grn.py   (158 ln)   Quantum (Ising) + NG-RC → GRN comparison
│
├── requirements.txt                   Pinned dependency list with explanations
├── setup.py                           setuptools install script (v0.2.0)
├── README.md                          User-facing documentation
└── MANIFEST.md                        This file
```

---

## Module Reference

### `classical/utils.py` — 324 lines

Data utilities and linear algebra helpers. All stateless functions.

| Function | Purpose |
|---|---|
| `preprocess(adata, n_hvg, log_transform)` | scanpy preprocessing: normalize → log1p → HVG → scale |
| `validate_expression_matrix(X)` | Check shape, finiteness; raises ValueError on failure |
| `order_by_pseudotime(X, pseudotime)` | Sort cell rows by pseudotime; returns X_sorted, pt_sorted, sort_idx |
| `compute_gene_velocity(X, pseudotime, window)` | Finite-difference dX/dt with moving-average smoothing |
| `ridge_regression(H, y, lambda_reg)` | Closed-form ridge: W = (H^T H + λI)^{-1} H^T y |
| `normalize_grn(GRN, method)` | Scale GRN to [0,1] via minmax / zscore / log |
| `threshold_grn(GRN, threshold)` | Binarize: 1 where GRN > threshold, else 0 |
| `compute_sensitivity(GRN)` | In-degree (normalized column sums) per gene |
| `compute_activity(GRN)` | Out-degree (normalized row sums) per gene |
| `get_top_regulators(GRN, gene_names, k)` | Top-k regulators per target gene |

---

### `classical/reservoir.py` — 277 lines

Builds fixed random weights and runs Echo State Network dynamics.

| Function | Purpose |
|---|---|
| `build_reservoir(n_reservoir, n_genes, spectral_radius, input_scale, density, random_state)` | Create W_res (sparse, scaled to spectral radius) and W_in (dense random). Returns weight dict. |
| `compute_reservoir_states(X, weights, leak_rate, washout)` | Run ESN dynamics over pseudotime-sorted cells. Returns **H matrix** (n_cells × n_reservoir). |
| `build_knn_graph(X, n_neighbors)` | Build Gaussian-weighted k-NN graph from expression data |
| `compute_graph_reservoir_states(X, weights, adjacency, leak_rate, washout, n_graph_iter)` | Reservoir states with iterative k-NN neighbor smoothing |

**Core dynamics equation** (inside `compute_reservoir_states`):
```
h(t+1) = (1 - α) · h(t)  +  α · tanh( W_res · h(t)  +  W_in · x(t) )
```
where α = `leak_rate`, x(t) = gene expression of pseudotime-ordered cell t.

---

### `classical/grn.py` — 252 lines

GRN inference by regression and back-projection. The only "training" in reservoir computing.

| Function | Purpose |
|---|---|
| `fit_readout_standard(H, X, lambda_reg)` | Regress X on H directly. No temporal structure. |
| `fit_readout_causal(H, X, pseudotime, lambda_reg)` | Regress X[t] on H[t-1]. Granger-causal, time-lagged. |
| `fit_readout_velocity(H, X, pseudotime, lambda_reg)` | Regress dX/dt on H midpoints. RNA-velocity style. |
| `infer_grn(W_in, W_out, normalize)` | GRN = \|W_in.T @ W_out\|. Back-projection through reservoir. |
| `infer_grn_from_dynamics(W_in, A, U, normalize)` | GRN from latent dynamics matrix A (landscape output). |
| `predict_expression(H, W_out)` | Predict gene expression: X_pred = H @ W_out |
| `compute_prediction_mse(X_true, H, W_out)` | Mean squared prediction error |

**GRN back-projection formula:**
```
GRN[i, j] = | Σ_k  W_in[k, i] · W_out[k, j] |
           = | (W_in.T @ W_out)[i, j] |
```
Gene i's influence on gene j = overlap between how i activates the reservoir (W_in) and how the reservoir predicts j (W_out).

---

### `classical/landscape.py` — 403 lines

Full Waddington attractor landscape pipeline. Each step is a separate function.

| Function | Purpose |
|---|---|
| `compress_reservoir_states(H, svd_rank, random_state)` | Randomized SVD: H ≈ U · diag(S) · Vt. Returns H_latent, U, S, Vt. |
| `compute_gene_velocity_landscape(X_sorted, pt_sorted)` | Finite-difference dX/dt (no smoothing, for landscape use) |
| `estimate_latent_dynamics(H_latent, pt_sorted, lambda_reg)` | Fit A in dH/dt ≈ H @ A via ridge regression on finite differences |
| `detect_attractors(X_sorted, velocity, n_attractors, random_state)` | Find low-velocity cells (bottom 10%), cluster with K-means → attractor centers |
| `compute_energy_landscape(H_latent, attractors)` | Energy = −log(KDE density). Low energy = attractor/valley. |
| `compute_fate_probabilities(H_latent, attractors)` | Sigmoid distance to each attractor, row-normalized to sum=1 |
| `run_landscape_pipeline(X, H, pseudotime, ...)` | Convenience wrapper: runs all 7 steps, returns full results dict |

**Full pipeline output dict** (from `run_landscape_pipeline`):

| Key | Shape | Meaning |
|---|---|---|
| `X_sorted` | (n_cells, n_genes) | Expression, pseudotime-sorted |
| `H_sorted` | (n_cells, n_reservoir) | Reservoir states, pseudotime-sorted |
| `pt_sorted` | (n_cells,) | Sorted pseudotime values |
| `sort_idx` | (n_cells,) | Indices used for sorting |
| `H_latent` | (n_cells, svd_rank) | SVD-compressed reservoir states |
| `U` | (n_cells, svd_rank) | Left singular vectors |
| `S` | (svd_rank,) | Singular values |
| `Vt` | (svd_rank, n_reservoir) | Right singular vectors |
| `velocity` | (n_cells, n_genes) | dX/dt per cell per gene |
| `A` | (svd_rank, svd_rank) | Latent dynamics matrix |
| `attractors` | (n_attractors, n_genes) | Attractor gene expression profiles |
| `cells_per_attractor` | list of arrays | Cell indices in each attractor basin |
| `energy` | (n_cells,) | Energy per cell (low = valley) |
| `fate_probs` | (n_cells, n_attractors) | Fate commitment probabilities |

---

### `classical/plotting.py` — 321 lines

All visualization functions. Accept numpy arrays directly — no class instances.

| Function | Purpose |
|---|---|
| `plot_grn_heatmap(GRN, gene_names, top_n, ax, cmap, figsize)` | Heatmap of top-n genes by out-degree |
| `plot_energy_landscape(energy, embedding, pseudotime, ax, figsize, cmap)` | 2D or 3D energy landscape scatter; auto-detects embedding shape |
| `plot_fate_probabilities(fate_probs, embedding, attractor_names, figsize)` | One subplot per attractor showing fate probability |
| `plot_attractor_genes(attractor_gene_programs, gene_names, top_n, ...)` | Heatmap of top variable genes per attractor |
| `plot_reservoir_dynamics(H, pseudotime, n_dims, figsize)` | PCA of H matrix, 2D or 3D, colored by pseudotime |

---

### `quantum/quantum_reservoir.py` — 870 lines

Full quantum reservoir computing implementation — classically simulated and NG-RC.
Organized in 9 sequential parts (A through I).

| Part | Contents |
|---|---|
| **A** — Primitives | `pauli_x/y/z()`, `identity()`, `kron_op_on_qubit()`, `kron_op_on_two_qubits()` |
| **B** — Hamiltonians | `build_ising_hamiltonian()`, `build_rrg_hamiltonian()`, `build_random_regular_graph()` |
| **C** — Unitary | `build_unitary(H, dt)` → U = exp(−iHΔt) |
| **D** — Encoding | `encode_input_rotation(psi, x_t, n_qubits, encoding_strength)` |
| **E** — Readout | `measure_observables(psi, n_qubits)`, `measure_density_observables(rho, n_qubits)` |
| **F** — Build weights | `build_quantum_reservoir(n_qubits, J, h, dt, ...)` → weight dict with H and U |
| **G** — Dynamics | `compute_quantum_reservoir_states(X, q_weights, ...)` → H_q matrix |
| **H** — NG-RC | `build_ngrc_features(X, delay_steps, poly_degree)`, `compute_ngrc_reservoir_states()` |
| **I** — Pipeline | `run_quantum_grn_pipeline(X, pseudotime, ...)` → full results dict |

**Quantum dynamics equation** (inside `compute_quantum_reservoir_states`):
```
|ψ(t+1)⟩ = U(Δt) · R_input(x_t) · |ψ(t)⟩
```
where U(Δt) = exp(−iHΔt) and R_input applies Y-rotations parameterized by x_t.

**NG-RC feature vector** (inside `build_ngrc_features`):
```
o_k  = [ x(t),  x(t−1), ...,  x(t−delay) ]     (time-delay embedding)
φ_k  = [ o_k ;  o_k ⊗ o_k ]                     (linear + all quadratic products)
```

---

## Feature Implementation Checklist

### Classical Reservoir (Echo State Network)
- [x] Sparse W_res with exact spectral radius scaling
- [x] Dense W_in, input_scale controlled
- [x] Leak-rate memory control (α)
- [x] Washout period to forget initial state
- [x] k-NN graph adjacency construction
- [x] Graph-smoothed reservoir states (n_graph_iter iterations)

### GRN Inference
- [x] Standard mode (no time structure)
- [x] Causal mode (time-lagged, Granger-style)
- [x] Velocity mode (dX/dt regression)
- [x] Back-projection: GRN = |W_in.T @ W_out|
- [x] Dynamics-based GRN from latent A matrix
- [x] Min-max / z-score / log normalization
- [x] Threshold sparsification
- [x] Top-k regulator ranking per gene
- [x] Prediction MSE evaluation

### Attractor Landscape
- [x] Randomized SVD compression (svd_rank controllable)
- [x] Finite-difference velocity computation
- [x] Ridge-regression latent dynamics (A matrix)
- [x] Low-velocity attractor detection (10th percentile)
- [x] K-means clustering of attractor candidates
- [x] KDE energy landscape with distance fallback
- [x] Sigmoid fate probabilities, row-normalized
- [x] Full pipeline convenience function

### Quantum Reservoir (Simulated)
- [x] Pauli matrices X, Y, Z and tensor product utilities
- [x] Transverse-field Ising Hamiltonian (1D chain, periodic)
- [x] Random Regular Graph Hamiltonian (from Paper 2)
- [x] Random regular graph adjacency generator
- [x] Unitary U = exp(−iHΔt) via scipy.linalg.expm
- [x] Input encoding via Y-rotation gates
- [x] Single-qubit observables ⟨Z_i⟩
- [x] Two-qubit correlators ⟨Z_i Z_j⟩
- [x] Density matrix observable measurement
- [x] Full quantum H_q matrix (n_cells × n_obs)

### Next-Generation RC (NG-RC)
- [x] Time-delay embedding o_k
- [x] Linear features (degree 1)
- [x] Quadratic features o_k ⊗ o_k (degree 2)
- [x] Cubic features (degree 3, for small gene sets)
- [x] Drop-in replacement for classical H matrix

### Utilities & Infrastructure
- [x] scanpy/AnnData compatibility
- [x] Reproducible random seeds throughout
- [x] Inline comments explaining every equation
- [x] No classes — all plain functions
- [x] setup.py with classical + quantum extras
- [x] requirements.txt with annotated dependencies

---

## Running the Examples

```bash
cd code/examples

# Classical pipeline: ESN → GRN on synthetic data
python 01_classical_grn.py

# Classical pipeline: ESN → attractor landscape + fate probabilities
python 02_classical_landscape.py

# Quantum pipeline: Ising reservoir vs NG-RC (both modes), GRN comparison
python 03_quantum_grn.py
```

Each script sets `sys.path` to find `classical/` and `quantum/` automatically
— no installation needed.

---

## Dependency Map

| Library | Used in | Reason |
|---|---|---|
| `numpy` | all modules | Array operations, matrix multiply |
| `scipy.sparse` | classical/reservoir.py | Sparse W_res (memory efficient) |
| `scipy.sparse.linalg.eigsh` | classical/reservoir.py | Spectral radius estimation |
| `scipy.linalg.expm` | quantum/quantum_reservoir.py | U = exp(−iHΔt) |
| `scipy.stats.gaussian_kde` | classical/landscape.py | Energy landscape density |
| `sklearn.neighbors.NearestNeighbors` | classical/reservoir.py | k-NN graph |
| `sklearn.cluster.KMeans` | classical/landscape.py | Attractor clustering |
| `sklearn.utils.extmath.randomized_svd` | classical/landscape.py | SVD compression |
| `sklearn.decomposition.PCA` | classical/plotting.py | Reservoir visualization |
| `scanpy` | classical/utils.py | Single-cell preprocessing |
| `anndata` | classical/utils.py | AnnData input format |
| `matplotlib` | classical/plotting.py | All base plots |
| `seaborn` | classical/plotting.py | Heatmaps |

---

## File Statistics

| File | Lines | Functions | Category |
|---|---|---|---|
| classical/utils.py | 324 | 10 | Utility |
| classical/reservoir.py | 277 | 4 | Core (classical) |
| classical/grn.py | 252 | 7 | Core (classical) |
| classical/landscape.py | 403 | 7 | Core (classical) |
| classical/plotting.py | 321 | 7 | Visualization |
| quantum/quantum_reservoir.py | 870 | 22 | Core (quantum) |
| examples/01_classical_grn.py | 111 | 2 | Example |
| examples/02_classical_landscape.py | 99 | 2 | Example |
| examples/03_quantum_grn.py | 158 | 2 | Example |
| **Total** | **2,815** | **63** | |

---

## Design Principles

1. **No classes** — every module is a collection of plain functions. Data flows explicitly through function arguments and return values. No hidden state.
2. **Sparse-first** — W_res is a `scipy.sparse.csr_matrix`. The k-NN adjacency matrix is also sparse. This keeps memory usage manageable for datasets with thousands of cells.
3. **Exact spectral radius** — W_res is rescaled so its largest eigenvalue (estimated via `eigsh` on W_res^T W_res) exactly equals the requested `spectral_radius`.
4. **Sequential clarity** — the landscape pipeline is split into 7 named functions that can each be called independently, inspected, and substituted.
5. **Quantum-classical symmetry** — `build_quantum_reservoir` mirrors `build_reservoir`; `compute_quantum_reservoir_states` mirrors `compute_reservoir_states`; both produce an H matrix usable with the same `fit_readout_*` and `infer_grn` functions.
6. **Reproducibility** — every random operation accepts a `random_state` seed.
