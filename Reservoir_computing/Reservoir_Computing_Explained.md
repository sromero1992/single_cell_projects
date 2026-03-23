# Reservoir Computing — Full Explanation
### Classical ESN · Quantum RC · scRNA-seq Application · Python Code Walkthrough
**scReservoir v0.2.0 — Selim Romero, TAMU**

---

## Table of Contents

1. Classical Reservoir Computing (Echo State Networks)
2. Pseudotime — How It Is Computed and Why It Matters
3. The H Matrix — The Bridge Between Cells and the Reservoir
4. Quantum Reservoir Computing
5. Next-Generation (NG-RC) Quantum-Inspired Features
6. Python Code — Every Function Explained
   - 6.1 `classical/utils.py`
   - 6.2 `classical/reservoir.py`
   - 6.3 `classical/grn.py`
   - 6.4 `classical/landscape.py`
   - 6.5 `classical/plotting.py`
   - 6.6 `quantum/quantum_reservoir.py`
7. The Full Pipeline in One Picture

---

## Part 1 — Classical Reservoir Computing

### 1.1 The Problem It Solves

Recurrent Neural Networks (RNNs) are the standard tool for sequential or time-series data. But training them with backpropagation through time (BPTT) is famously hard: gradients either vanish or explode as they travel backward through many time steps, making it nearly impossible to learn long-range dependencies.

**Reservoir Computing** removes this problem entirely by refusing to train the recurrent connections at all. The radical insight: pick random recurrent weights, fix them forever, and only train a simple linear layer on top. This sounds like it should fail — but for a huge class of problems it works just as well as, or better than, fully trained RNNs. The reason is that a large enough random recurrent network is already a powerful nonlinear feature extractor; you don't need to learn how to mix information, just how to read it out.

The two most influential variants are **Echo State Networks** (ESNs, Jaeger 2001) and **Liquid State Machines** (LSMs, Maass et al. 2002). scReservoir uses ESNs.

---

### 1.2 The Three-Layer Architecture

```
Input x(t) ──[W_in, fixed]──▶  Reservoir h(t) ──[W_out, TRAINED]──▶  Output y(t)
                                    ▲      │
                                    └──────┘  W_res (fixed, sparse, random)
```

**W_in** — Input weight matrix. Shape: `(n_reservoir × n_genes)`. Fixed random, drawn from Normal(0, input_scale). Projects gene expression into the high-dimensional reservoir space. Never modified.

**W_res** — Reservoir (recurrent) weight matrix. Shape: `(n_reservoir × n_reservoir)`. Fixed, sparse random. Only ~1% of entries are non-zero (`density=0.01`). Controls how the reservoir mixes information across neurons over time. Never modified.

**W_out** — Readout weight matrix. Shape: `(n_reservoir × n_genes)`. **The only thing trained.** Learned by ridge regression in one closed-form computation — no iterative gradient descent.

---

### 1.3 The Core Update Equation

At every time step `t` (each pseudotime-ordered cell), the reservoir state vector `h(t)` — of length `n_reservoir` — is updated:

```
h(t+1) = (1 − α) · h(t)  +  α · tanh( W_res · h(t)  +  W_in · x(t) )
```

Breaking this down term by term:

- **`x(t)`** is the gene expression vector of the cell at pseudotime step t. Shape: `(n_genes,)`.
- **`W_in · x(t)`** projects the current cell's expression into reservoir space. This is the "drive" from the input.
- **`W_res · h(t)`** is the recurrent contribution — the reservoir "remembers" its previous state and mixes it in. This creates temporal memory.
- **`tanh(...)`** squashes the combined signal into (−1, +1), creating the nonlinear transformation. Without this, the reservoir would just be a linear system and lose much of its power.
- **`α`** (leak rate) controls the balance between memory and update:
  - α close to 0: the state barely changes each step — long memory, slow response.
  - α = 1: the state is completely re-computed each step — no carry-over from the past.
  - Default: α = 0.3 (moderate memory).
- **`(1 − α) · h(t)`** is the "leaked" old state that carries over.

---

### 1.4 The Echo State Property and Spectral Radius

The **spectral radius** ρ of W_res is its largest absolute eigenvalue. It is the single most important hyperparameter.

The **Echo State Property** states: regardless of what initial state h(0) you start with, the effect of that starting condition must eventually "echo away" and be forgotten, so that the reservoir state is determined solely by the history of inputs. This is what gives the reservoir its memory of the input, not memory of its own arbitrary initial condition.

For the Echo State Property to hold, the spectral radius must generally be **< 1**. Why? Because if ρ(W_res) < 1, the recurrent dynamics are contractive — the reservoir state shrinks toward zero in the absence of input, guaranteeing that the initial condition is forgotten.

- **ρ < 1**: Stable, fading memory. The reservoir "echoes" past inputs for a while, then forgets them.
- **ρ ≈ 1**: Near the "edge of chaos." Maximum memory capacity and richest nonlinear dynamics, but the system becomes sensitive to noise.
- **ρ > 1**: Potentially unstable — reservoir states can grow without bound.

In code (`classical/reservoir.py`, `build_reservoir`):
1. Generate a random sparse matrix W_res.
2. Estimate its spectral radius ρ by computing the largest eigenvalue of W_res^T W_res via `eigsh`, then taking its square root.
3. Rescale: `W_res ← W_res × (spectral_radius / ρ)`.

The default is `spectral_radius = 0.9`.

---

### 1.5 Training — Only Ridge Regression

After running all cells through the reservoir to produce H, the only training step is:

```
Minimize:  ‖H · W_out − X‖²  +  λ · ‖W_out‖²
```

This has the closed-form solution:

```
W_out = (H^T H + λI)^{−1} H^T X
```

This is computed once, analytically, in a fraction of a second — no iterative optimization, no learning rate, no epochs, no backpropagation. λ (lambda_reg) prevents overfitting by penalizing large weights.

---

## Part 2 — Pseudotime: What It Is, How It Is Computed, and Why It Matters

### 2.1 The Single-Cell Data Problem

When you perform single-cell RNA sequencing on a tissue, you capture a **snapshot** of many cells at the same moment. If you are studying development — say, how a stem cell becomes a neuron — you cannot watch the same cell over time (the sequencing process destroys the cell). Instead you capture thousands of cells that are at different stages of the developmental process at the time of capture.

The challenge is: these cells have no time labels. You have a cloud of gene expression vectors with no indication of which ones are "early" and which are "late" in development.

**Pseudotime** is a computational solution: assign each cell a scalar value between 0 and 1 that estimates how far along the developmental trajectory it is. A cell with pseudotime ≈ 0 is likely near the stem/progenitor state; a cell with pseudotime ≈ 1 is likely near a terminal, differentiated state.

### 2.2 The Biological Intuition

The reason pseudotime works is that development is a **continuous process**. Even though you have a snapshot, cells that are at similar developmental stages will have similar gene expression profiles. Cells at nearby pseudotime values will cluster together in the gene expression space. The developmental trajectory traces a path (or branching paths) through this high-dimensional space.

Finding pseudotime is equivalent to finding the ordering of points along this path.

### 2.3 Step 1 — Dimensionality Reduction

Single-cell data has thousands of genes (dimensions). Working directly in this space is computationally intractable and dominated by noise. The first step is to reduce dimensionality while preserving the structure of the data.

The standard pipeline:
1. **Normalize** — divide each cell's counts by the total counts, then multiply by 10,000. This removes the effect of how deeply each cell was sequenced (sequencing depth normalization).
2. **Log-transform** — apply log(1 + count). Gene expression data is highly skewed; this compression makes the data more symmetric and amenable to analysis.
3. **Select HVGs** — keep only the top 2,000 Highly Variable Genes (genes whose expression varies the most across cells). These carry the most information about cell state differences.
4. **PCA** — Principal Component Analysis reduces from 2,000 dimensions to ~50 principal components. These components capture the major directions of variation, which typically correspond to biological processes (cell type, developmental stage, cell cycle).

### 2.4 Step 2 — Build a k-NN Graph

From the PCA coordinates, construct a **k-Nearest Neighbor graph**: connect each cell to its k (typically 15–30) most similar neighbors in PCA space. This graph captures the local topology of the data — which cells are developmentally adjacent.

Critically, the k-NN graph is computed in PCA space (not the full gene space), which is more robust to noise.

### 2.5 Step 3 — Graph Partitioning into Clusters (Leiden/Louvain)

Before computing pseudotime, it is useful to identify the major cell populations (clusters). The Leiden or Louvain algorithm finds communities in the k-NN graph — groups of cells that are more densely connected to each other than to the rest. These communities correspond to biologically distinct cell types or states.

This step is optional for pseudotime but helps define the root cell (where pseudotime = 0).

### 2.6 Step 4 — UMAP for Visualization

**UMAP** (Uniform Manifold Approximation and Projection) is a nonlinear dimensionality reduction that takes the k-NN graph and computes a 2D embedding where nearby cells in graph space are placed nearby in 2D. This is primarily used for visualization — you can see the trajectory as a curve or tree in the UMAP plot.

UMAP itself does not assign pseudotime — it only reveals the topology visually.

### 2.7 Step 5 — Pseudotime via Diffusion Pseudotime (DPT)

The most widely used pseudotime algorithm in single-cell biology is **Diffusion Pseudotime (DPT)** by Haghverdi et al. (2016), implemented as `sc.tl.dpt()` in scanpy.

**The core idea**: imagine a random walk on the k-NN graph. Starting from a "root cell" (the most progenitor-like cell, manually selected or automatically detected), how many steps does it take for the random walk to reach each other cell?

Cells that are reachable in few steps are early (low pseudotime). Cells that require many steps are late (high pseudotime).

Formally, DPT uses **diffusion maps** — a spectral method that computes eigenvectors of the normalized transition matrix of the k-NN graph. The diffusion distances in this spectral space are used to order cells.

The algorithm:
1. Construct a normalized transition matrix T from the k-NN graph (T[i,j] = probability of walking from cell i to cell j in one step).
2. Compute the top eigenvectors of T (via eigendecomposition).
3. Embed cells in the space of these eigenvectors — this is the "diffusion map embedding."
4. The pseudotime of each cell is the diffusion distance from the root cell in this embedding.

**Other pseudotime methods** include:
- **Monocle 3** (DDRTree / principal graphs): fits a "principal tree" through the data and measures arc-length along the tree.
- **Slingshot**: fits smooth curves (principal curves) through clusters in dimensionality-reduced space and measures arc-length.
- **PAGA** (Partition-based Graph Abstraction): computes connectivity between Leiden clusters to build a coarse-grained trajectory graph, then interpolates positions within each cluster.

### 2.8 Step 6 — Choosing the Root Cell

For DPT, you must specify a root cell — the cell most likely to be at the start of the developmental trajectory. Strategies for finding it:

- **Biological knowledge**: if you know which cluster corresponds to stem cells or progenitors, pick a cell from that cluster.
- **RNA velocity**: the direction of mRNA splicing provides a temporal arrow. The cells where all arrows point "away" are likely roots.
- **Diffusion component extrema**: the cell with the smallest or largest value along the first diffusion component is often near the trajectory endpoint.

In scanpy:
```python
adata.uns['iroot'] = np.argmin(adata.obsm['X_diffmap'][:, 0])   # example root selection
sc.tl.dpt(adata)
pseudotime = adata.obs['dpt_pseudotime'].values
```

### 2.9 Why Pseudotime Is Central to scReservoir

Reservoir computing is a **temporal framework** — it was designed for time-series data where inputs arrive in a meaningful temporal order. The reservoir accumulates memory of past inputs and uses that memory to make predictions.

Single-cell snapshots have no natural time order. Pseudotime creates that order. By sorting cells by pseudotime before feeding them to the reservoir, you are treating the developmental trajectory as a time series: cell at pseudotime 0 is time step 1, cell at pseudotime 0.001 is time step 2, and so on.

The reservoir then builds up a running memory of the entire trajectory: the reservoir state at cell t encodes not just cell t's gene expression, but a nonlinear, fading memory of every earlier cell in pseudotime. This is what makes causal GRN inference possible — earlier gene expression changes are "remembered" in the reservoir state when later expression is predicted.

**A key limitation**: pseudotime is an estimate, not a measurement. It assumes the developmental trajectory is approximately tree-like and that the k-NN graph topology reflects developmental proximity. For complex, branching, or circular trajectories, pseudotime can be misleading, and the causal GRN inferred by scReservoir will reflect the estimated ordering, not a true temporal one.

### 2.10 Pseudotime in the Code

In `classical/utils.py`:

```python
X_sorted, pt_sorted, sort_idx = order_by_pseudotime(X, pseudotime)
```

This is literally just `np.argsort(pseudotime)` applied to the rows of X. Simple. The complexity is entirely in how pseudotime itself was computed upstream (by scanpy/DPT, before scReservoir is called).

Then in `classical/reservoir.py`, `compute_reservoir_states(X_sorted, weights)` processes the sorted cells one by one in the for loop, building up the H matrix row by row.

---

## Part 3 — The H Matrix

### 3.1 What It Is

After running all n_cells pseudotime-sorted cells through the reservoir dynamics, you collect the reservoir state vector h(t) after each step into a matrix H:

```
H =  [ h(1) ]      ← reservoir state after seeing cell 1  (pseudotime ≈ 0)
     [ h(2) ]      ← reservoir state after seeing cell 2
     [ h(3) ]      ← reservoir state after seeing cell 3
     [  ...  ]
     [ h(T) ]      ← reservoir state after seeing cell T  (pseudotime ≈ 1)
```

**H has shape `(n_cells × n_reservoir)`**. In the default configuration this is `(n_cells × 500)`.

### 3.2 Why H Is Not Just a Nonlinear Transform of X

The crucial point that makes reservoir computing powerful is that **H[t] is not just a function of X[t]**. Because of the recurrent term `W_res · h(t)` in the update equation, H[t] depends on all previous reservoir states, which in turn depend on all previous inputs. So H[t] is a nonlinear function of the entire sequence X[0], X[1], ..., X[t].

In other words: H[50] encodes something about what gene expression looked like at pseudotime steps 0, 1, 2, ..., 49 — in addition to what it looks like at step 50. The older information is there but "faded" (because of the spectral radius < 1 and the leak rate < 1).

This is the reservoir's "memory" — and it is what allows causal GRN inference.

### 3.3 The Washout Period

The initial reservoir state is h(0) = 0 (all zeros). This is arbitrary and has nothing to do with the biology. During the first ~100 steps (the washout period), the reservoir state is dominated by this arbitrary starting point rather than by the input data.

To handle this, the code runs `washout` extra steps before starting to record H. During washout, cells are cycled through repeatedly but none of the states are recorded. By the time washout ends, the effect of h(0) = 0 has been "echoed away" and the reservoir state is determined by the input data.

```python
for t in range(n_cells + washout):
    if t < washout:
        x_t = X[t % n_cells]    # cycle through data
    else:
        x_t = X[t - washout]    # actual cell in order
    # ... update h ...
    if t >= washout:
        H[t - washout] = h      # only record after washout
```

---

## Part 4 — Quantum Reservoir Computing

### 4.1 What Changes

In classical reservoir computing the "reservoir" is a network of artificial neurons. In **Quantum Reservoir Computing (QRC)**, the reservoir is a **quantum many-body system** — a collection of quantum bits (qubits) evolving under a fixed Hamiltonian.

The conceptual structure is identical:

| Classical ESN | Quantum Reservoir |
|---|---|
| W_res (random, fixed recurrent weights) | Hamiltonian H (fixed, determines dynamics) |
| W_in (projects input into reservoir) | Encoding gates R_input(x) (rotates qubits) |
| h(t+1) = (1-α)h(t) + α·tanh(W_res·h(t) + W_in·x(t)) | \|ψ(t+1)⟩ = U(Δt) · R_input(x(t)) · \|ψ(t)⟩ |
| h(t) ∈ ℝ^{n_reservoir} | \|ψ(t)⟩ ∈ ℂ^{2^n_qubits} |
| Readout: h(t) directly | Readout: ⟨Z_i⟩, ⟨Z_i Z_j⟩ measured |
| W_out trained by ridge regression | W_out trained by ridge regression (identical) |

The key difference is in what the "reservoir state" is and how rich it can be. An n-qubit quantum state lives in a 2^n-dimensional Hilbert space. With 10 qubits you have a 1024-dimensional state space — more than most classical reservoirs — achieved with just 10 physical components.

### 4.2 The Hamiltonian

The Hamiltonian H describes the energy and dynamics of the quantum system. scReservoir implements two:

**Transverse-field Ising model** (Paper 1 — Sornsaeng et al.):
```
H = −J · Σ_{i} Z_i Z_{i+1}  +  h · Σ_{i} X_i
```
- Z_i, X_i are Pauli operators on qubit i (2×2 matrices acting on one qubit, tensored with identities on all others)
- J = Ising coupling (favors aligned neighboring spins)
- h = transverse field (tilts spins toward the X axis, creating quantum superpositions)
- The competition between J and h drives a quantum phase transition near h/J ≈ 1 — the "critical point"
- **Near the critical point, the system has maximum correlations and memory — optimal for reservoir computing** (the quantum analog of spectral radius ≈ 1 at the edge of chaos)

**Random Regular Graph Hamiltonian** (Paper 2 — Ivaki et al.):
```
H = Σ_{(i,j)∈E} [J_z · Z_i Z_j + J_x · X_i X_j]  +  Σ_i [(h_z)·Z_i + (h_x + δ_i)·X_i]
```
- E = edge set of a random regular graph (each qubit has exactly k neighbors)
- δ_i = random local fields (disorder): drives the system from integrable to chaotic
- Key finding: **optimal performance is at the boundary between chaotic and localized phases** (edge of chaos, again)
- Graph degree k controls information propagation speed: higher k = faster spread, but too dense loses locality

### 4.3 Time Evolution — The Unitary

The quantum reservoir evolves for one time step Δt via the unitary operator:

```
U(Δt) = exp(−i H Δt)
```

This is the quantum analog of W_res. Just as W_res drives h(t) → h(t+1), U drives |ψ(t)⟩ → |ψ(t+1)⟩. Unlike W_res, U is exactly **unitary**: it preserves the norm of the quantum state. This is guaranteed because H is Hermitian (H = H^†), which makes exp(−iHΔt) unitary.

In code: `scipy.linalg.expm(-1j * H * dt)`.

### 4.4 Input Encoding

Each cell's gene expression x(t) is injected into the quantum state by applying Y-rotation gates to the qubits:

```
R_y(θ_i) = [[cos(θ_i/2), −sin(θ_i/2)],
             [sin(θ_i/2),  cos(θ_i/2)]]
```

where θ_i = encoding_strength × x_t[bin_i]. Since n_genes >> n_qubits, genes are averaged into n_qubit bins first. The full rotation operator is R = R_y(θ_0) ⊗ R_y(θ_1) ⊗ ... ⊗ R_y(θ_{n-1}).

After encoding: |ψ'⟩ = R_input(x_t) · |ψ(t)⟩. Then after evolution: |ψ(t+1)⟩ = U · |ψ'⟩.

### 4.5 Readout — Measuring Observables

The quantum state |ψ(t+1)⟩ cannot be read out directly (it is a complex vector of length 2^n_qubits). Instead, you measure **expectation values of quantum observables**. scReservoir measures:

- **Single-qubit Z**: ⟨Z_i⟩ = ⟨ψ|Z_i|ψ⟩ for i = 0..n_qubits−1. This is the average magnetization along Z for qubit i.
- **Two-qubit ZZ correlators**: ⟨Z_i Z_j⟩ for all pairs i < j. This captures quantum correlations (entanglement) between qubits.

Total observables per time step: n_qubits + n_qubits(n_qubits−1)/2. For n_qubits=6: 6 + 15 = 21.

These real-valued numbers form one row of the H_q matrix — the quantum analog of h(t).

### 4.6 Why Quantum Is Different (Not Just Bigger)

A classical reservoir with 21 neurons would be trivially weak. But 21 quantum observables from 6 qubits carry information about a state in a 64-dimensional Hilbert space. The key advantages:

**Entanglement**: the two-qubit correlators ⟨Z_i Z_j⟩ carry quantum correlations that cannot be decomposed into products of single-qubit properties. These non-classical correlations provide a form of nonlinear mixing that is provably hard to simulate classically.

**Exponential state space**: even though you only read out 21 numbers, the underlying quantum state explored 64 dimensions to produce them. Each observable is a weighted average over this exponential space.

**Inherent temporal dynamics**: quantum unitary evolution is exact (no floating-point decay) and provides perfect memory of past inputs as long as the system remains coherent.

---

## Part 5 — Next-Generation (NG-RC) Quantum-Inspired Features

From Paper 1 (Sornsaeng et al.):

**What if you could get quantum-like nonlinear mixing without a quantum computer?**

The NG-RC insight: the key thing a quantum reservoir does is compute high-order products of past input values through entanglement. You can approximate this classically by explicitly constructing polynomial features of a time-delay embedding:

```
o_k = [ x(t),  x(t−1),  ...,  x(t−m+1) ]        (time-delay embedding, length m × n_genes)
φ_k = [ o_k ;  o_k ⊗ o_k ]                        (linear + all pairwise products)
```

The outer product `o_k ⊗ o_k` computes ALL pairwise products of current and past gene expression values. For m=2 and n_genes=40: o_k has length 80, and o_k ⊗ o_k has length 6400. The full feature vector φ_k has 6480 dimensions, encoding every possible pairwise interaction between gene expression at time t and gene expression at time t−1.

This is what quantum entanglement achieves physically — correlating current and past values through joint quantum states. The NG-RC approach computes it explicitly and classically.

**Advantages**: no quantum simulation, scales to any dataset size, interpretable features.
**Disadvantage**: feature dimension grows as O(n_genes² × m²), which can become large for many genes.

In the code (`quantum/quantum_reservoir.py`, Part H): `compute_ngrc_reservoir_states(X, delay_steps=1, poly_degree=2)`. The output H_ng matrix is used with the exact same `fit_readout_*` and `infer_grn()` functions as classical and quantum H matrices.

---

## Part 6 — Python Code: Every Function Explained

### 6.1 `classical/utils.py`

**`preprocess(adata, n_hvg, log_transform)`**

Standard single-cell preprocessing pipeline using scanpy. Four steps in sequence:
1. `sc.pp.normalize_total(target_sum=1e4)` — each cell's counts sum to 10,000 (removes sequencing depth differences)
2. `sc.pp.log1p()` — log(1 + counts) compress the skewed distribution
3. `sc.pp.highly_variable_genes(n_top_genes=n_hvg)` — keep the top 2,000 most variable genes (discard housekeeping genes)
4. `sc.pp.scale(max_value=10)` — standardize each gene to zero mean / unit variance; clip outliers at 10

Returns a dense float32 numpy array of shape (n_cells, n_hvg).

---

**`validate_expression_matrix(X)`**

Three guard checks before running the pipeline:
- `X.ndim != 2` → raise ValueError ("must be 2D")
- `X.shape[0] == 0 or X.shape[1] == 0` → raise ValueError ("empty dimensions")
- `not np.isfinite(X).all()` → raise ValueError ("NaN or Inf")

Returns True if all pass.

---

**`order_by_pseudotime(X, pseudotime)`**

Sorts the rows of X according to ascending pseudotime:
```python
sort_idx = np.argsort(pseudotime)
return X[sort_idx], pseudotime[sort_idx], sort_idx
```
Returns: X_sorted, pt_sorted, sort_idx. The sort_idx is returned so you can apply the same ordering to other arrays (cell labels, metadata, etc.).

---

**`compute_gene_velocity(X, pseudotime, window)`**

Estimates dX/dt for each cell using finite differences:
```
velocity[t] ≈ (X[t] − X[t−1]) / (pseudotime[t] − pseudotime[t−1])
```
The first cell is set equal to the second (no predecessor). Then applies a moving-average smoother (window=5) to reduce noise from the finite-difference approximation: `velocity[t] ← mean(velocity[t−w/2 : t+w/2])`.

Returns shape (n_cells, n_genes). This is the gene-level velocity — how fast each gene is being expressed at each pseudotime point.

---

**`ridge_regression(H, y, lambda_reg)`**

Closed-form ridge regression. Solves:
```
min_W  ‖H @ W − y‖²  +  λ ‖W‖²
```
Via: `W = (H^T H + λI)^{−1} H^T y`

Implemented step by step:
1. `HTH = H.T @ H` — Gram matrix
2. `HTy = H.T @ y` — cross-correlation
3. `HTH_reg = HTH + lambda_reg * I` — Tikhonov regularization
4. `W = np.linalg.solve(HTH_reg, HTy)` — solve the linear system (more stable than explicit inverse)

---

**`normalize_grn(GRN, method)`**, **`threshold_grn(GRN, threshold)`**

`normalize_grn`: scale to [0,1] via minmax (default), zscore, or log1p+scale. Always clips final values to [0,1].

`threshold_grn`: returns a binary (0/1) matrix. Entry [i,j] = 1 if GRN[i,j] > threshold.

---

**`compute_sensitivity(GRN)`**, **`compute_activity(GRN)`**, **`get_top_regulators(GRN, gene_names, k)`**

`compute_sensitivity`: normalized column sums of GRN. GRN[:,j] = all influences flowing INTO gene j. High sensitivity = heavily regulated gene (many inputs).

`compute_activity`: normalized row sums of GRN. GRN[i,:] = all influences flowing OUT of gene i. High activity = master regulator / transcription factor.

`get_top_regulators`: for each target gene j, sorts all genes by GRN[:,j] (column j) and returns the top-k. Returns a dict: `{gene_j: [(regulator_1, score), (regulator_2, score), ...]}`.

---

### 6.2 `classical/reservoir.py`

**`build_reservoir(n_reservoir, n_genes, spectral_radius, input_scale, density, random_state)`**

The entire reservoir is constructed here once and stored in a plain dict. Never called again.

For W_res:
1. `scipy.sparse.random(n_reservoir, n_reservoir, density=density)` — sparse matrix with `density` fraction of non-zero entries
2. `eigsh(W_res^T W_res, k=1)` → largest eigenvalue of W_res^T W_res → `ρ = sqrt(eigenvalue)`
3. `W_res *= spectral_radius / ρ` — rescale to exact target spectral radius

For W_in:
- `rng.randn(n_reservoir, n_genes) * input_scale` — dense Gaussian random matrix

Returns a dict: `{'W_res': ..., 'W_in': ..., 'n_reservoir': ..., 'n_genes': ..., 'spectral_radius': ..., ...}`.

---

**`compute_reservoir_states(X, weights, leak_rate, washout)`**

Runs the full ESN dynamics over n_cells pseudotime-ordered cells and produces H.

The loop runs `n_cells + washout` total iterations:
- First `washout` iterations: cycle through data (`X[t % n_cells]`), update h, but do NOT record. This "primes" the reservoir.
- Remaining `n_cells` iterations: process cells in order, record h after each step.

Update at each step:
```python
drive = W_res @ h + W_in @ x_t
h = (1 - leak_rate) * h  +  leak_rate * np.tanh(drive)
```

Returns H of shape (n_cells, n_reservoir). Each row is the reservoir's internal state snapshot after processing one cell.

---

**`build_knn_graph(X, n_neighbors)`**

Builds a weighted k-NN adjacency matrix from gene expression:
1. `NearestNeighbors(n_neighbors).fit(X).kneighbors(X)` — find k nearest neighbors for each cell
2. Edge weight = `exp(−distance²)` — Gaussian kernel, so nearby cells get strong connections
3. Returns a `scipy.sparse.csr_matrix` of shape (n_cells, n_cells)

---

**`compute_graph_reservoir_states(X, weights, adjacency, leak_rate, washout, n_graph_iter)`**

Two-step process:
1. Run standard `compute_reservoir_states(X, weights)` to get initial H
2. Iteratively smooth H across the k-NN graph `n_graph_iter` times:
   ```python
   A_norm = adjacency / row_sums    # row-normalize
   H = np.tanh(A_norm @ H)          # propagate + nonlinear activation
   ```

The effect: each cell's reservoir state is blended toward its transcriptionally similar neighbors. Cells in the same developmental state end up with more similar H rows. This is a form of graph-based regularization.

---

### 6.3 `classical/grn.py`

**`fit_readout_standard(H, X, lambda_reg)`**

Simplest mode. Calls `ridge_regression(H, X, lambda_reg)` directly.

Assumption: every row of H and X are paired (same cell ordering). No temporal structure is enforced.

Returns W_out of shape (n_reservoir, n_genes).

---

**`fit_readout_causal(H, X, pseudotime, lambda_reg)`**

Granger-causal mode. Sorts both H and X by pseudotime, then creates a time lag:
```python
H_lagged = H_sorted[:-1]    # "past" reservoir states (cells 0 to T-2)
X_target = X_sorted[1:]     # "future" gene expression (cells 1 to T-1)
W_out = ridge_regression(H_lagged, X_target, lambda_reg)
```

The logic: if gene A causes gene B, then A's signal should appear in the reservoir state (H) before B's expression rises in X. By regressing X[t] on H[t−1], we force the model to predict the future from the past, making causal relationships learnable.

This is analogous to Granger causality in time-series econometrics.

---

**`fit_readout_velocity(H, X, pseudotime, lambda_reg)`**

RNA-velocity-inspired mode. Instead of predicting expression levels, predicts the rate of change:
```python
velocity = dX / dt           # shape (n_cells-1, n_genes)
H_mid = (H[:-1] + H[1:]) / 2     # midpoint reservoir states
W_out = ridge_regression(H_mid, velocity, lambda_reg)
```

Using dX/dt as the target is more sensitive to directional dynamics — a gene that is rising at time t will have positive velocity, even if its absolute level is still low. This can reveal regulatory relationships that are invisible in snapshot-based approaches.

---

**`infer_grn(W_in, W_out, normalize)`**

The key inference step. Computes:
```
GRN = |W_in.T @ W_out|         shape: (n_genes × n_genes)
```

The reasoning:
- `W_in[k, i]` = how strongly gene i activates reservoir neuron k
- `W_out[k, j]` = how strongly reservoir neuron k predicts gene j
- `(W_in.T @ W_out)[i, j] = Σ_k W_in[k,i] · W_out[k,j]` = the total influence of gene i on gene j mediated through all reservoir neurons

Taking the absolute value ensures all influences are non-negative. Then normalize_grn scales to [0,1].

---

**`infer_grn_from_dynamics(W_in, A, U, normalize)`**

Alternative GRN using the dynamics matrix A from the landscape module:
```python
A_expanded = U @ A @ U.T            # project A from latent (svd_rank) to full reservoir space
GRN = |W_in.T @ A_expanded @ W_in|  # project through input/output weights
```

A encodes how the latent reservoir state evolves over pseudotime. By projecting through W_in (which maps genes → reservoir), you get a gene-gene influence matrix that reflects the temporal dynamics of the trajectory, not just the instantaneous regression.

---

**`predict_expression(H, W_out)`** and **`compute_prediction_mse(X_true, H, W_out)`**

`predict_expression`: `H @ W_out` — straightforward matrix multiply.

`compute_prediction_mse`: `mean((X_true - H @ W_out)²)` — average squared prediction error over all cells and all genes.

---

### 6.4 `classical/landscape.py`

**`compress_reservoir_states(H, svd_rank, random_state)`**

Performs randomized SVD:
```
H ≈ U · diag(S) · Vt
```
- n_components = min(svd_rank, n_cells, n_reservoir)
- Uses `sklearn.utils.extmath.randomized_svd` — a fast approximate algorithm that does not compute the full SVD (much faster for large matrices)
- H_latent = U @ diag(S), shape (n_cells, svd_rank)

Think of this as PCA of the reservoir states: H_latent is the "scores matrix" (where each cell lives in the low-dimensional latent space), and Vt is the "loadings matrix" (which reservoir neurons define each latent dimension).

---

**`compute_gene_velocity_landscape(X_sorted, pt_sorted)`**

Simpler version of velocity (no smoothing). Finite differences:
```
velocity[t] = (X[t] − X[t−1]) / (pt[t] − pt[t−1])
```
First cell set = second cell. Returns shape (n_cells, n_genes).

---

**`estimate_latent_dynamics(H_latent, pt_sorted, lambda_reg)`**

Fits the linear system: `dH_latent/dt ≈ H_latent @ A`

Uses finite differences of H_latent and ridge regression:
```python
dt  = np.diff(pt_sorted)
dH  = np.diff(H_latent, axis=0)
H_prev = H_latent[:-1]
# Ridge: A^T = (H_prev^T H_prev + λI)^{-1} H_prev^T (dH/dt)
A_T = solve(H_prev^T @ H_prev + λI, H_prev^T @ (dH / dt))
A = A_T.T
```

A is a (svd_rank × svd_rank) matrix — the "dynamics matrix." The eigenvalues of A describe the stability of the trajectory:
- Eigenvalues with negative real part → converging dynamics (approaching an attractor)
- Eigenvalues with positive real part → diverging (unstable)

---

**`detect_attractors(X_sorted, velocity, n_attractors, random_state)`**

Finds stable states — cells that have stopped changing:
1. Compute per-cell velocity magnitude: `||velocity[t]||` (L2 norm over genes)
2. Flag cells in the bottom 10th percentile: `vel_magnitude < np.percentile(vel_magnitude, 10)`
3. Cluster these low-velocity cells with `KMeans(n_clusters=n_attractors)`
4. Return: cluster centers (n_attractors × n_genes) as `attractors`, and per-cluster cell indices as `cells_per_attractor`

The intuition: a cell that has reached a terminal differentiated state has "settled" — its gene expression is no longer changing. Low velocity in pseudotime = developmental stability = attractor.

---

**`compute_energy_landscape(H_latent, attractors)`**

Assigns an energy value (landscape height) to every cell:
```python
kde    = gaussian_kde(H_latent.T)     # fit KDE to all cell positions in latent space
energy = -log(kde(H_latent.T) + ε)   # negative log-density
```

KDE estimates the probability density p(h) at each cell's latent position. Dense regions (where many cells cluster) have high p(h) → low energy. Sparse regions (transient states) have low p(h) → high energy.

The result: valleys in the landscape (low energy, high density) correspond to attractors / terminal cell states. Hills correspond to transitional unstable states.

Fallback: if the KDE covariance matrix is singular (degenerate data), uses summed Euclidean distance to attractor centers instead.

---

**`compute_fate_probabilities(H_latent, attractors)`**

For each cell and each attractor:
1. Compute Euclidean distance from cell's latent position to attractor center: `d = ||H_latent[t] − attractor_k||`
2. Convert to probability via sigmoid: `p = exp(−d) / (1 + exp(−d))`
3. Row-normalize so probabilities across all attractors sum to 1

Returns fate_probs of shape (n_cells, n_attractors). Cell t's probability of ending up at attractor k.

---

**`run_landscape_pipeline(X, H, pseudotime, ...)`**

Convenience function that runs all 7 landscape steps in sequence and returns every intermediate result in one dict. Useful for interactive analysis where you want to inspect each step.

---

### 6.5 `classical/plotting.py`

All plotting functions accept numpy arrays directly (no class instances). They all accept an optional `ax` parameter — if provided, they draw onto it; otherwise they create a new figure.

**`plot_grn_heatmap(GRN, gene_names, top_n, ax, cmap, figsize)`**

Selects top_n genes by out-degree (row sum of GRN), extracts the top_n × top_n submatrix, and draws a seaborn heatmap. Rows = regulators, columns = targets.

**`plot_energy_landscape(energy, embedding, pseudotime, ax, figsize, cmap)`**

Dispatches to `_plot_energy_2d` or `_plot_energy_3d` based on `embedding.shape[1]`. Scatter plot colored by energy value. If pseudotime is provided, draws a dashed line connecting cells in pseudotime order (the developmental trajectory).

**`plot_fate_probabilities(fate_probs, embedding, attractor_names, figsize)`**

Creates one subplot per attractor. Each subplot is a scatter plot colored by the probability of reaching that attractor. All subplots share the same [0,1] color scale.

**`plot_attractor_genes(attractor_gene_programs, gene_names, top_n, ...)`**

Selects top_n genes by variance across attractors (genes that are most different between attractors). Draws a heatmap with attractors on x-axis and genes on y-axis.

**`plot_reservoir_dynamics(H, pseudotime, n_dims, figsize)`**

Runs PCA on H to get 2D or 3D coordinates. Scatter plot colored by pseudotime. Draws the pseudotime trajectory as a dashed line. Axis labels show explained variance ratio for each PC.

---

### 6.6 `quantum/quantum_reservoir.py`

**Part A — Quantum Primitives**

`pauli_x()`, `pauli_y()`, `pauli_z()`, `identity(d)`: return the standard 2×2 (or d×d) complex numpy arrays.

`kron_op_on_qubit(op, qubit_idx, n_qubits)`: embeds a single-qubit operator on qubit `qubit_idx` inside an n_qubit system by tensoring with identity on all other qubits. For example, Z on qubit 1 of a 3-qubit system: I ⊗ Z ⊗ I.

`kron_op_on_two_qubits(op_i, op_j, qubit_i, qubit_j, n_qubits)`: same for two-qubit operators. Used to construct ZZ interaction terms.

---

**Part B — Hamiltonians**

`build_ising_hamiltonian(n_qubits, J, h, periodic)`:

Constructs the full (2^n × 2^n) Hamiltonian matrix:
```
H = −J · Σ_i ZZ_{i,i+1}  +  h · Σ_i X_i
```
Loops over all bonds (i, i+1) and all sites i, adding the corresponding tensor-product operators. `periodic=True` adds the wrap-around bond Z_{n-1} Z_0.

`build_rrg_hamiltonian(n_qubits, adjacency, J_z, J_x, h_z, h_x, disorder_strength, random_state)`:

Same construction but uses the adjacency matrix instead of a 1D chain. Also adds random local fields `δ_i ∈ [−Δ, Δ]` that break translational symmetry and drive the integrable-to-chaotic transition.

`build_random_regular_graph(n_qubits, degree, random_state)`:

Constructs a k-regular graph (every node has exactly k neighbors) by first creating a ring with degree/2 connections per side, then randomly rewiring edges while preserving the degree sequence.

---

**Part C — Unitary**

`build_unitary(H, dt)`: computes `U = scipy.linalg.expm(-1j * H * dt)`. This is the matrix exponential of a Hermitian matrix, giving a unitary matrix. For n_qubits=6, H is 64×64, and expm is fast. For n_qubits=10, H is 1024×1024 — still manageable on CPU.

---

**Part D — Input Encoding**

`encode_input_rotation(psi, x_t, n_qubits, encoding_strength)`:

Maps n_genes gene expression values to n_qubits rotation angles (by averaging into bins), builds the rotation operator R = ⊗_i R_y(θ_i), and applies it to the quantum state: `psi_encoded = R @ psi`.

---

**Part E — Readout**

`measure_observables(psi, n_qubits)`:

For each qubit i: `⟨Z_i⟩ = Re(psi^† @ Z_i_full @ psi)` where `Z_i_full` is the full n_qubit Z operator on qubit i.
For each pair (i,j): `⟨Z_i Z_j⟩ = Re(psi^† @ ZiZj_full @ psi)`.
Returns array of length n_qubits + n_qubits(n_qubits−1)/2.

`measure_density_observables(rho, n_qubits)`: same but using `np.trace(rho @ observable)` for density matrices.

---

**Part F — Build Quantum Reservoir**

`build_quantum_reservoir(n_qubits, J, h, dt, hamiltonian_type, ...)`:

Analogous to `build_reservoir`. Calls the appropriate Hamiltonian builder, then `build_unitary`. Returns a dict: `{'H': ..., 'U': ..., 'n_qubits': ..., 'n_obs': ..., 'dt': ..., ...}`.

---

**Part G — Quantum Dynamics (produces H_q)**

`compute_quantum_reservoir_states(X, q_weights, encoding_strength, washout)`:

The quantum analog of `compute_reservoir_states`. Same loop structure:
1. Washout phase: cycle through data, update |ψ⟩, do not record.
2. Recording phase: for each cell t: (a) encode x_t via R_input, (b) apply U, (c) renormalize |ψ⟩, (d) measure observables → row of H_q.

Returns H_q of shape (n_cells, n_obs). H_q[t] = [⟨Z_i(t)⟩, ⟨Z_i Z_j(t)⟩] after processing cells 0 through t.

---

**Part H — NG-RC Features**

`build_ngrc_features(X, delay_steps, poly_degree)`:

For each valid cell t (from delay_steps to n_cells−1):
1. Build time-delay embedding: `o = concat([x(t), x(t−1), ..., x(t−delay_steps)])`, length (delay_steps+1) × n_genes
2. Append linear features: φ = [o]
3. If poly_degree ≥ 2: append `np.outer(o, o).flatten()` — all pairwise products
4. If poly_degree ≥ 3: append triple products (very large — use sparingly)

The first delay_steps cells are dropped (no prior context). The matching expression matrix for regression is `X[delay_steps:]`.

---

**Part I — Full Quantum Pipeline**

`run_quantum_grn_pipeline(X, pseudotime, ...)`:

End-to-end convenience function. Sorts X by pseudotime, builds the quantum reservoir or NG-RC features, fits readout weights, infers GRN, and returns everything in a dict. The `mode` argument ('standard', 'causal', 'velocity') is implemented inline since it needs to align the quantum H_q with the correct rows of X.

---

## Part 7 — The Full Pipeline in One Picture

```
Raw scRNA-seq data (AnnData)
         │
         ▼
  utils.preprocess(adata)
    → normalize + log1p + HVG selection + scale
    → X: (n_cells × n_HVG)
         │
         ▼
  Compute pseudotime externally (scanpy DPT):
    sc.pp.neighbors → k-NN graph
    sc.tl.diffmap   → diffusion map embedding
    sc.tl.dpt       → pseudotime per cell
    → pseudotime: (n_cells,), values ∈ [0,1]
         │
         ▼
  utils.order_by_pseudotime(X, pseudotime)
    → X_sorted: cells ordered from progenitor → terminal
    → pt_sorted: sorted pseudotime values
         │
         ├─────────────────────────────────────────────────┐
         │  CLASSICAL PATH                                  │  QUANTUM PATH
         ▼                                                  ▼
  reservoir.build_reservoir(n_reservoir=500)    quantum_reservoir.build_quantum_reservoir(n_qubits=6)
    → weights = {W_res (sparse), W_in (dense)}    → q_weights = {H (Hamiltonian), U (unitary)}
         │                                                  │
         ▼                                                  ▼
  reservoir.compute_reservoir_states(X_sorted)  quantum_reservoir.compute_quantum_reservoir_states(X_sorted)
  OR: reservoir.compute_graph_reservoir_states  OR: quantum_reservoir.compute_ngrc_reservoir_states
    → H: (n_cells × 500)                          → H_q: (n_cells × 21)
         │                                                  │
         └──────────────────┬──────────────────────────────┘
                            │  Both produce an H matrix used identically downstream
                            ▼
  grn.fit_readout_causal(H, X_sorted, pt_sorted)
    → regress X[t] on H[t−1]  (Granger causality)
    → W_out: (n_reservoir × n_genes)
         │
         ▼
  grn.infer_grn(weights['W_in'], W_out)
    → GRN = |W_in.T @ W_out|        (back-projection)
    → GRN: (n_genes × n_genes)     normalized to [0,1]
         │
         ├──────────────────────────────────────────────────
         │  Also run landscape analysis:
         ▼
  landscape.run_landscape_pipeline(X_sorted, H, pt_sorted)
    → Step 1: compress_reservoir_states → H_latent (SVD)
    → Step 2: compute_gene_velocity_landscape → velocity
    → Step 3: estimate_latent_dynamics → A matrix
    → Step 4: detect_attractors → attractor centers
    → Step 5: compute_energy_landscape → energy per cell
    → Step 6: compute_fate_probabilities → fate_probs
    → Returns full results dict
         │
         ▼
  Visualization:
    plotting.plot_grn_heatmap(GRN, gene_names)
    plotting.plot_energy_landscape(energy, UMAP_embedding)
    plotting.plot_fate_probabilities(fate_probs, UMAP_embedding)
    plotting.plot_attractor_genes(attractor_programs, gene_names)
    plotting.plot_reservoir_dynamics(H, pseudotime)
```

---

## Key Equations Reference

| Name | Formula | Where Used |
|---|---|---|
| ESN update | h(t+1) = (1−α)h(t) + α·tanh(W_res·h(t) + W_in·x(t)) | reservoir.compute_reservoir_states |
| Ridge regression | W = (H^T H + λI)^{-1} H^T y | utils.ridge_regression |
| GRN back-projection | GRN = \|W_in.T @ W_out\| | grn.infer_grn |
| SVD compression | H ≈ U·diag(S)·Vt | landscape.compress_reservoir_states |
| Latent dynamics | dH/dt ≈ H @ A | landscape.estimate_latent_dynamics |
| KDE energy | energy = −log(p_KDE(h)) | landscape.compute_energy_landscape |
| Fate probability | p = exp(−d) / (1 + exp(−d)), normalized | landscape.compute_fate_probabilities |
| Quantum evolution | \|ψ(t+1)⟩ = U·R(x_t)·\|ψ(t)⟩ | quantum_reservoir.compute_quantum_reservoir_states |
| Quantum unitary | U(Δt) = exp(−iHΔt) | quantum_reservoir.build_unitary |
| Ising Hamiltonian | H = −J·ΣZZ + h·ΣX | quantum_reservoir.build_ising_hamiltonian |
| NG-RC features | φ = [o; o⊗o] where o = [x(t), x(t−1)] | quantum_reservoir.build_ngrc_features |
| Observable | ⟨Z_i⟩ = Re(ψ^\dagger Z_i ψ) | quantum_reservoir.measure_observables |

---

*Document last updated: scReservoir v0.2.0*
