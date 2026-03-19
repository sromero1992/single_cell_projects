# scReservoir — MATLAB Toolbox

Reservoir computing framework for single-cell RNA-seq analysis.
Infers gene regulatory networks (GRNs), reconstructs Waddington energy landscapes, and detects cell fate attractors.

---

## Requirements

- MATLAB R2020b or later
- Toolboxes: Statistics and Machine Learning, Signal Processing (optional)
- No additional toolboxes required for core functionality

---

## Installation

```matlab
addpath(genpath('/path/to/scReservoir/matlab'));
```

Or add the path permanently via MATLAB's `pathtool`.

---

## Quick Start

```matlab
% 1. Preprocess
[X, genes, ~] = scReservoir_preprocess(X_raw, geneNames, 'nHVG', 2000);

% 2. Initialize reservoir
res = scReservoir_init(size(X,2), 'N_res', 500, 'rho', 0.9);

% 3. Infer GRN (causal mode)
[GRN, topRegs] = scReservoir_GRN(res, X, pseudotime, genes, 'mode', 'causal', 'k', 10);

% 4. Display top regulators
disp(topRegs(1))

% 5. Full landscape pipeline
result = scReservoir_landscape(X, pseudotime, genes, 'nAttractors', 5);
```

---

## Function Reference

| Function | Description |
|---|---|
| `scReservoir_preprocess` | Normalize, log-transform, select HVGs |
| `scReservoir_init` | Initialize random reservoir (W_res, W_in) |
| `scReservoir_run` | Drive reservoir with expression input, return state matrix H |
| `scReservoir_GRN` | GRN inference: static, causal, or velocity mode |
| `scReservoir_graphGRN` | Graph-based reservoir (kNN cell graph propagation) |
| `scReservoir_scalable` | Atlas-scale: sparse reservoir + randomized SVD |
| `scReservoir_landscape` | Full pipeline: dynamics + attractors + landscape + GRN |
| `scReservoir_plot_GRN` | Visualize GRN as heatmap and network graph |
| `scReservoir_plot_landscape` | Visualize Waddington landscape and fate probabilities |

---

## GRN Inference Modes

```matlab
% Static: regress expression on reservoir states
[GRN, ~] = scReservoir_GRN(res, X, pt, genes, 'mode', 'static');

% Causal: lagged reservoir enforces temporal directionality
[GRN, ~] = scReservoir_GRN(res, X, pt, genes, 'mode', 'causal');

% Velocity: regress dx/dt — models continuous-time gene regulation
[GRN, ~] = scReservoir_GRN(res, X, pt, genes, 'mode', 'velocity');

% Graph: propagate across kNN cell graph (handles branching)
[GRN, ~, H] = scReservoir_graphGRN(res, X, pt, genes, 'kNN', 15);

% Ensemble averaging (reduces noise)
[GRN, ~] = scReservoir_GRN(res, X, pt, genes, 'mode', 'causal', 'ensemble', 10);
```

---

## Scalable Mode (>50k cells)

```matlab
[GRN, topRegs, Hlatent] = scReservoir_scalable(X, pseudotime, genes, ...
    'N_res',    800, ...
    'rankSVD',  80, ...
    'density',  0.01, ...
    'mode',     'velocity');
```

Recommended settings by dataset size:

| Dataset | N_res | SVD rank | Density |
|---|---|---|---|
| 10k cells | 400 | 40 | 0.05 |
| 100k cells | 600 | 60 | 0.02 |
| 1M cells | 800–1200 | 80–120 | 0.01 |

---

## Waddington Landscape

```matlab
result = scReservoir_landscape(X, pseudotime, genes, ...
    'N_res',       600, ...
    'rankSVD',     50, ...
    'nAttractors', 5, ...
    'plot',        true);

% Access outputs
result.GRN               % nGenes x nGenes regulatory matrix
result.attractorCells    % indices of attractor cells
result.fateProbabilities % nCells x nAttractors
result.energy            % Waddington potential per cell
result.pca_scores        % 2D PCA of latent states
result.topRegulators     % top regulators per gene
```

---

## Mathematical Framework

**Reservoir dynamics (leaky integrator):**
```
h(t) = (1-α)*h(t-1) + α*tanh(W_res*h(t-1) + W_in*x(t))
```

**GRN back-projection:**
```
influence(i,g) = |W_in' * W_out^g|_i
```

**Velocity mode (continuous-time GRN):**
```
dx/dt ≈ H * W_out  →  GRN(i,g) = |W_in' * W_out^g|_i
```

**Latent dynamics:**
```
dh/dt = A*h  →  GRN = |W_in' * A * W_in|
```

**Waddington potential:**
```
U(h) ≈ -log P(h)   (KDE in PCA space)
```

---

## Examples

| Script | Description |
|---|---|
| `examples/example01_basic_GRN.m` | Static and causal GRN from scratch |
| `examples/example02_graph_and_velocity.m` | Graph reservoir + velocity mode |
| `examples/example03_attractor_landscape.m` | Full landscape pipeline with fate probabilities |

---

## Citation

If you use scReservoir in your research, please cite:

> [Authors]. scReservoir: Reservoir Dynamical Systems for Gene Regulatory Network Inference and Cell Fate Attractor Reconstruction from Single-Cell Transcriptomics. *Bioinformatics*, [year].

Also cite the reservoir computing paper that inspired this work:

> Prediction performance of random reservoirs with input-driven dynamics. *Chaos: An Interdisciplinary Journal of Nonlinear Science*, 36(3), 033117 (2026).
