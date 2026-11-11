"""
Attractor landscape functions for scReservoir.

Models the Waddington developmental landscape through a sequential pipeline
of plain functions.  Each function does exactly one thing and returns its
result explicitly — no side effects, no shared state.

Typical usage (full pipeline)
------------------------------
H_latent, U, S, Vt = compress_reservoir_states(H, svd_rank=50)
velocity            = compute_gene_velocity(X_sorted, pt_sorted)
A                   = estimate_latent_dynamics(H_latent, pt_sorted)
attractors, cells_per_attractor = detect_attractors(X_sorted, velocity, n_attractors=5)
energy              = compute_energy_landscape(H_latent, attractors)
fate_probs          = compute_fate_probabilities(H_latent, attractors)
"""

import numpy as np
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
from sklearn.utils.extmath import randomized_svd


# ---------------------------------------------------------------------------
# Step 1 — Compress reservoir states  (SVD / PCA-like reduction)
# ---------------------------------------------------------------------------

def compress_reservoir_states(H, svd_rank=50, random_state=42):
    """
    Compress reservoir state matrix H via randomized SVD.

    The 500-neuron reservoir state space is high-dimensional and noisy.
    SVD finds the directions of maximum variance and projects H into
    a much smaller space (svd_rank dimensions), analogous to PCA.

    Decomposition:
        H  ≈  U  @  diag(S)  @  Vt

    The compressed (latent) representation used downstream is:
        H_latent  =  U  @  diag(S)        shape (n_cells, svd_rank)

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix.  Must already be sorted by pseudotime.
    svd_rank : int, default=50
        Number of singular vectors to keep.  Higher = more detail, more noise.
    random_state : int, default=42
        Seed for the randomized SVD algorithm.

    Returns
    -------
    H_latent : np.ndarray of shape (n_cells, svd_rank)
        Cells in compressed latent space.
    U : np.ndarray of shape (n_cells, svd_rank)
        Left singular vectors (cell scores).
    S : np.ndarray of shape (svd_rank,)
        Singular values (importance of each latent dimension).
    Vt : np.ndarray of shape (svd_rank, n_reservoir)
        Right singular vectors (how reservoir neurons contribute).
    """
    n_components = min(svd_rank, H.shape[0], H.shape[1])
    U, S, Vt = randomized_svd(H, n_components=n_components, random_state=random_state)
    H_latent = U @ np.diag(S)
    return H_latent, U, S, Vt


# ---------------------------------------------------------------------------
# Step 2 — Gene velocity  (finite differences along pseudotime)
# ---------------------------------------------------------------------------

def compute_gene_velocity_landscape(X_sorted, pt_sorted):
    """
    Compute gene velocity via finite differences along sorted pseudotime.

    Estimates dX/dt for each cell as the change in expression divided by
    the change in pseudotime between adjacent cells.

    Note: This is a simpler version without moving-average smoothing,
    intended for the landscape pipeline.  For smoothed velocity, see
    compute_gene_velocity() in utils.py.

    Parameters
    ----------
    X_sorted : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix, sorted by ascending pseudotime.
    pt_sorted : np.ndarray of shape (n_cells,)
        Sorted pseudotime values.

    Returns
    -------
    velocity : np.ndarray of shape (n_cells, n_genes)
        Finite-difference velocity.  The first cell is given the same
        velocity as the second cell (no predecessor to difference against).
    """
    dt = np.diff(pt_sorted)                      # (n_cells-1,)
    dX = np.diff(X_sorted, axis=0)               # (n_cells-1, n_genes)

    velocity = np.zeros_like(X_sorted)
    velocity[1:] = dX / dt[:, None]
    velocity[0]  = velocity[1]                   # replicate second cell's velocity

    return velocity


# ---------------------------------------------------------------------------
# Step 3 — Latent dynamics  (fit a linear dynamical system to H_latent)
# ---------------------------------------------------------------------------

def estimate_latent_dynamics(H_latent, pt_sorted, lambda_reg=1e-3):
    """
    Fit a linear dynamics matrix A that describes how the latent state evolves.

    Model:  dH_latent / dt  ≈  H_latent  @  A

    A is fit by ridge regression on finite-difference estimates of the
    time derivative of H_latent.  A encodes the "flow" of the developmental
    trajectory in compressed reservoir space.

    The eigenvalues of A determine the stability:
      - Negative real parts  → stable attractor (converging dynamics)
      - Positive real parts  → unstable / diverging dynamics

    Parameters
    ----------
    H_latent : np.ndarray of shape (n_cells, svd_rank)
        Compressed reservoir states from compress_reservoir_states().
        Must be sorted by pseudotime.
    pt_sorted : np.ndarray of shape (n_cells,)
        Sorted pseudotime values.
    lambda_reg : float, default=1e-3
        Regularization strength for ridge regression.

    Returns
    -------
    A : np.ndarray of shape (svd_rank, svd_rank)
        Latent dynamics matrix.  Satisfies  dH/dt ≈ H @ A.
    """
    dt = np.diff(pt_sorted)                         # (n_cells-1,)
    dH = np.diff(H_latent, axis=0)                  # (n_cells-1, svd_rank)

    H_prev = H_latent[:-1]                          # states at time t

    # Weighted ridge regression:  minimize ||dH/dt - H_prev @ A||² + λ||A||²
    HTH  = H_prev.T @ H_prev
    HTdH = H_prev.T @ (dH / dt[:, None])            # weight each pair by 1/dt

    A_T = np.linalg.solve(
        HTH + lambda_reg * np.eye(HTH.shape[0]),
        HTdH
    )
    A = A_T.T

    return A


# ---------------------------------------------------------------------------
# Step 4 — Attractor detection  (low-velocity cells = stable states)
# ---------------------------------------------------------------------------

def detect_attractors(X_sorted, velocity, n_attractors=5, random_state=42):
    """
    Identify attractor states as regions of low gene-expression velocity.

    The intuition: a cell that has stopped changing (low dX/dt magnitude) has
    reached a stable state — an attractor in the Waddington landscape.
    Conversely, cells still changing rapidly are mid-transition.

    Algorithm:
      1. Compute per-cell velocity magnitude: ||velocity[t]||
      2. Flag cells in the bottom 10th percentile as attractor candidates.
      3. K-means cluster these candidates into n_attractors groups.
      4. Return cluster centers as the attractor gene expression profiles.

    Parameters
    ----------
    X_sorted : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (pseudotime-sorted).
    velocity : np.ndarray of shape (n_cells, n_genes)
        Gene velocities from compute_gene_velocity_landscape().
    n_attractors : int, default=5
        Number of distinct attractor states (e.g., cell fates) to look for.
    random_state : int, default=42
        Seed for K-means.

    Returns
    -------
    attractors : np.ndarray of shape (n_attractors, n_genes)
        Gene expression centroid of each attractor.
    cells_per_attractor : list of np.ndarray
        cells_per_attractor[i] = indices of cells belonging to attractor i.
    """
    vel_magnitude = np.linalg.norm(velocity, axis=1)    # (n_cells,)

    # Cells with low velocity are potential attractor members
    threshold    = np.percentile(vel_magnitude, 10)
    attractor_mask  = vel_magnitude < threshold
    attractor_cells = np.where(attractor_mask)[0]

    if len(attractor_cells) > n_attractors:
        # Cluster the low-velocity cells to find distinct attractor states
        kmeans = KMeans(n_clusters=n_attractors, random_state=random_state, n_init=10)
        labels = kmeans.fit_predict(X_sorted[attractor_cells])

        attractors = kmeans.cluster_centers_                    # (n_attractors, n_genes)
        cells_per_attractor = [
            attractor_cells[labels == i] for i in range(n_attractors)
        ]
    else:
        # Too few low-velocity cells — use them all as one group
        attractors = X_sorted[attractor_cells]
        cells_per_attractor = [attractor_cells]

    return attractors, cells_per_attractor


# ---------------------------------------------------------------------------
# Step 5 — Energy landscape  (KDE-based; low energy = attractor valleys)
# ---------------------------------------------------------------------------

def compute_energy_landscape(H_latent, attractors=None):
    """
    Assign an "energy" value to every cell, analogous to Waddington landscape height.

    Energy is estimated as the negative log-density of the cell distribution
    in latent space (kernel density estimation):
        energy[t]  =  -log( p(H_latent[t]) )

    Dense regions (where many cells accumulate) → low energy (valleys / attractors).
    Sparse regions (transient, unstable states) → high energy (hills / saddle points).

    If KDE fails (singular covariance matrix), falls back to using the summed
    Euclidean distance from each cell to the detected attractor centers.

    Parameters
    ----------
    H_latent : np.ndarray of shape (n_cells, svd_rank)
        Compressed reservoir states from compress_reservoir_states().
    attractors : np.ndarray of shape (n_attractors, n_genes), optional
        Attractor gene expression profiles from detect_attractors().
        Used only as a fallback when KDE fails.

    Returns
    -------
    energy : np.ndarray of shape (n_cells,)
        Energy value per cell.  Lower = more stable / attractor-like.
    """
    if H_latent.shape[0] < 10:
        return np.zeros(H_latent.shape[0])

    try:
        kde    = gaussian_kde(H_latent.T)          # KDE fit to all cells
        energy = -np.log(kde(H_latent.T) + 1e-10)  # negative log-density
    except np.linalg.LinAlgError:
        # Fallback: use sum of distances to attractor centers
        energy = np.zeros(H_latent.shape[0])
        if attractors is not None:
            for attractor in attractors:
                dist = np.linalg.norm(H_latent - attractor[np.newaxis, :], axis=1)
                energy += dist

    return energy


# ---------------------------------------------------------------------------
# Step 6 — Fate probabilities  (distance-to-attractor → probability)
# ---------------------------------------------------------------------------

def compute_fate_probabilities(H_latent, attractors):
    """
    Compute the probability that each cell will commit to each attractor.

    For each attractor, a cell's proximity in latent space is converted to
    a probability using a sigmoid-like function (exp(-distance)), then
    probabilities across all attractors are normalized to sum to 1.

    Parameters
    ----------
    H_latent : np.ndarray of shape (n_cells, svd_rank)
        Compressed reservoir states from compress_reservoir_states().
    attractors : np.ndarray of shape (n_attractors, n_genes)
        Attractor gene expression profiles from detect_attractors().

    Returns
    -------
    fate_probs : np.ndarray of shape (n_cells, n_attractors)
        fate_probs[t, k] = probability that cell t will reach attractor k.
        Rows sum to 1.
    """
    n_cells      = H_latent.shape[0]
    n_attractors = len(attractors)

    fate_probs = np.zeros((n_cells, n_attractors))

    for k, attractor in enumerate(attractors):
        dist = np.linalg.norm(H_latent - attractor[np.newaxis, :], axis=1)
        # Sigmoid-like: closer cells get higher probability
        fate_probs[:, k] = np.exp(-dist) / (1.0 + np.exp(-dist))

    # Normalize across attractors so each row sums to 1
    row_sums = fate_probs.sum(axis=1, keepdims=True)
    fate_probs = fate_probs / (row_sums + 1e-10)

    return fate_probs


# ---------------------------------------------------------------------------
# Convenience: run the full landscape pipeline in one call
# ---------------------------------------------------------------------------

def run_landscape_pipeline(X, H, pseudotime, n_attractors=5, svd_rank=50, lambda_reg=1e-3):
    """
    Run the complete attractor landscape analysis in sequence.

    This is a convenience wrapper that calls all six steps in order and
    returns every intermediate result in a single dictionary.  You can also
    call each step independently if you want to inspect or modify the
    intermediate outputs.

    Steps performed:
      1. Sort X and H by pseudotime
      2. Compress H via SVD  → H_latent, U, S, Vt
      3. Compute gene velocity  → velocity
      4. Estimate latent dynamics  → A
      5. Detect attractors  → attractors, cells_per_attractor
      6. Compute energy landscape  → energy
      7. Compute fate probabilities  → fate_probs

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (unsorted OK — will be sorted here).
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix (unsorted OK — will be sorted here).
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime per cell.
    n_attractors : int, default=5
        Number of attractor states to identify.
    svd_rank : int, default=50
        SVD compression rank.
    lambda_reg : float, default=1e-3
        Regularization for dynamics estimation.

    Returns
    -------
    results : dict with keys:
        'X_sorted'            : (n_cells, n_genes)  expression, sorted by pseudotime
        'H_sorted'            : (n_cells, n_reservoir)  states, sorted by pseudotime
        'pt_sorted'           : (n_cells,)  sorted pseudotime
        'sort_idx'            : (n_cells,)  indices used for sorting
        'H_latent'            : (n_cells, svd_rank)  compressed states
        'U'                   : (n_cells, svd_rank)
        'S'                   : (svd_rank,)
        'Vt'                  : (svd_rank, n_reservoir)
        'velocity'            : (n_cells, n_genes)
        'A'                   : (svd_rank, svd_rank)  latent dynamics matrix
        'attractors'          : (n_attractors, n_genes)
        'cells_per_attractor' : list of arrays
        'energy'              : (n_cells,)
        'fate_probs'          : (n_cells, n_attractors)
    """
    # Step 1: sort by pseudotime
    sort_idx  = np.argsort(pseudotime)
    X_sorted  = X[sort_idx]
    H_sorted  = H[sort_idx]
    pt_sorted = pseudotime[sort_idx]

    # Step 2: compress reservoir states
    H_latent, U, S, Vt = compress_reservoir_states(H_sorted, svd_rank=svd_rank)

    # Step 3: compute gene velocity
    velocity = compute_gene_velocity_landscape(X_sorted, pt_sorted)

    # Step 4: estimate latent dynamics
    A = estimate_latent_dynamics(H_latent, pt_sorted, lambda_reg=lambda_reg)

    # Step 5: detect attractors
    attractors, cells_per_attractor = detect_attractors(
        X_sorted, velocity, n_attractors=n_attractors
    )

    # Step 6: energy landscape
    energy = compute_energy_landscape(H_latent, attractors)

    # Step 7: fate probabilities
    fate_probs = compute_fate_probabilities(H_latent, attractors)

    return {
        'X_sorted':            X_sorted,
        'H_sorted':            H_sorted,
        'pt_sorted':           pt_sorted,
        'sort_idx':            sort_idx,
        'H_latent':            H_latent,
        'U':                   U,
        'S':                   S,
        'Vt':                  Vt,
        'velocity':            velocity,
        'A':                   A,
        'attractors':          attractors,
        'cells_per_attractor': cells_per_attractor,
        'energy':              energy,
        'fate_probs':          fate_probs,
    }
