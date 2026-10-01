# What CellRank Did in Script 10 — A Plain-Language Walkthrough

This document explains the trajectory-inference method used in
`10_trajectory_cellrank.R`, step by step, including the pieces that were not
obvious: the CytoTRACE pseudotime, the smoothed `Ms` layer (the "moments"
approximation), the Markov transition matrix, and how terminal states and fate
probabilities are derived. It describes *this* pipeline's exact configuration,
not CellRank in the abstract.

---

## 1. The big idea

CellRank models differentiation as a **Markov chain over cells**. Every cell is
a state; the probability of "moving" from one cell to a neighbouring cell
encodes the direction of differentiation. Once you have that transition matrix,
a lot follows almost for free: where cells end up (terminal states), where they
start (initial states), the probability each cell reaches each endpoint (fate
probabilities), an ordering along the process (pseudotime), and a projected
"flow field" you can draw as arrows.

The crucial question is: **where does the direction come from?** Different
CellRank *kernels* answer this differently. RNA-velocity kernels use
spliced/unspliced dynamics. We do **not** have velocity data (standard Cell
Ranger output has no spliced/unspliced counts), so we use the
**CytoTRACEKernel**, which derives direction from a CytoTRACE-style potency
score instead. No velocity required.

---

## 2. CytoTRACE, and the pseudotime it produces

**CytoTRACE** (Gulati et al., *Science* 2020) is built on one empirical
observation: **the number of distinct genes a cell expresses tracks its
developmental potential.** Stem/progenitor cells keep many transcriptional
programs simultaneously accessible, so they express *more* genes; as a cell
commits and matures, expression concentrates into a focused program and the
gene count drops.

CellRank's CytoTRACEKernel reimplements this:

1. For each cell it counts expressed genes (the "gene-count signature").
2. It finds the genes whose expression correlates most with that signature.
3. Using a **smoothed** expression layer (see §3), it computes a per-cell
   `ct_score` — high = more stem-like, low = more differentiated.
4. It defines a pseudotime: **`ct_pseudotime = 1 − ct_score`**. So **low
   pseudotime = less differentiated (earlier)**, high = more differentiated.

This is the "pseudotime inference" you noticed: it is not Monocle/DPT diffusion
pseudotime — it is a potency-derived ordering. In our object these land in the
metadata as `cellrank_ct_score`, `cellrank_ct_pseudotime`,
`cellrank_ct_num_exp_genes`.

---

## 3. The `Ms` layer — the "moments" / smoothing approximation

Raw single-cell counts are extremely sparse (mostly zeros from dropout).
Correlating each gene against the gene-count signature on raw counts would be
dominated by that dropout noise. CytoTRACE therefore needs a **smoothed**
expression matrix.

In the RNA-velocity world (scvelo), smoothing is done by computing **moments**:
`Ms` = the first moment (kNN-average) of the *spliced* counts, and `Mu` = the
first moment of the *unspliced* counts. "First moment" just means: replace each
cell's value with the average over its nearest neighbours. This is the `Xmu`/
`Ms`/`Mu` machinery you were asking about.

Two things matter for us:

- The CytoTRACEKernel only needs the **smoothed expression** (`Ms`). It does
  **not** need `Mu`, because our directionality comes from gene counts, not from
  splicing dynamics.
- `scv.pp.moments()` refuses to run without spliced/unspliced layers (it prints
  *"Skipping moments, because un/spliced counts were not found"*). We don't have
  those layers.

So the script builds `Ms` **directly** from scanpy's neighbour graph, which is
exactly what a first moment is — a one-hop kNN average of the log-normalised
expression:

```
Ms = rownormalize(connectivities + I) @ X
```

where `connectivities` is the cell–cell kNN graph, `+ I` includes each cell
itself, `rownormalize` makes each row sum to 1 (so it's a weighted average), and
`X` is the log-normalised matrix. That single line is the "approximation of
`Xmu`" — the smoothed layer the kernel reads from. It is standard and correct;
it just sidesteps scvelo's velocity-only requirement.

---

## 4. The transition matrix

With `ct_pseudotime` in hand, the CytoTRACEKernel builds a **row-stochastic
transition matrix** `T`: for each cell, it looks at its kNN neighbours and
assigns higher transition probability to neighbours that lie **forward** along
the pseudotime gradient (more differentiated), and lower probability to those
behind it. The result is a directed random walk that, on average, flows from
stem-like toward differentiated cells.

This is what `ctk$compute_transition_matrix()` produces. Everything downstream
is analysis *of this matrix*.

---

## 5. Macrostates, terminal states, initial states (GPCCA)

A 30,000×30,000 transition matrix is too fine-grained to interpret directly, so
CellRank **coarse-grains** it into a handful of **macrostates** — metastable
groups of cells the random walk tends to linger in. The algorithm is **GPCCA**
(Generalized Perron Cluster Cluster Analysis), which uses a **real Schur
decomposition** of the transition matrix.

- **Schur spectrum (`schur_spectrum.png`)** — the eigenvalues of `T`.
  Eigenvalues near 1 correspond to slow, persistent (metastable) modes; the
  number of them *before the largest gap* is the natural number of macrostates.
  This is why the script plots the spectrum first: you set `N_STATES` from that
  gap rather than guessing.
- **Macrostates** (`g$compute_macrostates(n_states = N_STATES)`) — the coarse
  states themselves, named by their dominant `CellType` (that's what
  `CLUSTER_KEY` controls).
- **Terminal states** (`g$predict_terminal_states(n_states = N_TERMINAL_STATES)`)
  — the subset of macrostates identified as differentiation **endpoints**
  (absorbing/most-stable states). Initial states are the opposite end.

Because these are statistical constructs, the script includes a **potency
sanity check** (STEP 7b): a true terminal state should sit at LOW CytoTRACE2
potency (most differentiated). If one lands on high-potency cells, it's flagged
as suspicious.

---

## 6. Fate probabilities

Given the terminal states, CellRank computes, **for every cell, the probability
that a random walk started there ends in each terminal state**
(`g$compute_fate_probabilities()`). Mathematically these are **absorption
probabilities** of the Markov chain — solved as a linear system on the
transition matrix. They land in the object as `fate_<name>` columns, plus a
`fate_dominant` (most likely endpoint) and `fate_confidence` (its probability).

---

## 7. The projected flow (arrows / streamlines)

The transition matrix lives in high-dimensional space. To visualise it,
`ctk$plot_projection()` projects each cell's expected next-step displacement down
onto the **UMAP** (the same `umap_harmony` embedding carried over from Seurat),
giving a 2D vector per cell. Those vectors are rendered as **streamlines**
(`flow_streamlines.png`) and **grid arrows** (`flow_arrows.png`). Arrows point
toward increasing differentiation — they should flow *out* of the stem/TA
compartment toward the mature endpoints.

---

## 8. Exactly how this pipeline was configured

| Choice | Value | Why |
|---|---|---|
| Kernel | CytoTRACEKernel | no RNA velocity available |
| Smoothing (`Ms`) | manual kNN average of log-norm X | scvelo moments needs spliced/unspliced |
| Embedding for neighbours | reused Harmony (`X_harmony`, 20 dims) | batch-corrected, so batch ≠ trajectory |
| HVGs / PCs / neighbours | 2000 / 20 / 30 | standard graph resolution |
| Cells | downsampled to 30,000 | GPCCA densifies the matrix; memory-bounded |
| `N_STATES` / `N_TERMINAL_STATES` | set from `schur_spectrum.png` | data-driven, not guessed |
| Threads | 1 (serial) | avoids parallel workers duplicating big arrays |
| Cross-check | terminal states vs CytoTRACE2 potency | catch wrong endpoint calls |

---

## 9. What to remember when interpreting results

- **`ct_pseudotime` is potency-derived**, not diffusion pseudotime. Low = less
  differentiated.
- **Terminal states are statistical endpoints, not proven biology.** Confirm
  them with known markers and the potency cross-check before drawing
  conclusions.
- The **direction is inferred from gene counts (CytoTRACE)**, not from spliced/
  unspliced dynamics — so it reflects a potency gradient, which is appropriate
  for a clean differentiation hierarchy (stem/TA → colonocytes, goblet, etc.)
  but should not be over-read as mechanistic velocity.
- If the flow arrows point *into* the stem compartment, or a terminal state is
  flagged as high-potency, revisit `N_STATES`/`N_TERMINAL_STATES` (via the
  spectrum) or the embedding.

---

*Key references:* Lange et al., *CellRank* (Nat. Methods 2022) and CellRank 2
(Weiler et al. 2024); Gulati et al., *CytoTRACE* (Science 2020); Reuter et al.,
GPCCA (2018–2019); Bergen et al., *scVelo* moments (Nat. Biotech 2020).
