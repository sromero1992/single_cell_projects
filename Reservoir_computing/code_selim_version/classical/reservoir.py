"""
Reservoir computing functions for scReservoir.

Implements Echo State Networks (ESNs) as plain functions.
All state (weights, etc.) is stored in plain dicts and passed explicitly
— no classes, no hidden state, no side effects.

Typical usage
-------------
# 1. Build the fixed random weights once
weights = build_reservoir(n_reservoir=500, n_genes=X.shape[1])

# 2. Run the reservoir dynamics to get the H matrix
H = compute_reservoir_states(X_sorted, weights)

# 3. (optional) Graph-regularized version
adjacency = build_knn_graph(X, n_neighbors=10)
H = compute_graph_reservoir_states(X_sorted, weights, adjacency)
"""

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------------------
# Weight initialization
# ---------------------------------------------------------------------------

def build_reservoir(
    n_reservoir,
    n_genes,
    spectral_radius=0.9,
    input_scale=0.5,
    density=0.01,
    random_state=42
):
    """
    Build the two fixed (never-trained) weight matrices of an Echo State Network.

    The recurrent weight matrix W_res is sparse and randomly generated.
    Its largest eigenvalue (spectral radius) is scaled to the requested value.
    The input weight matrix W_in is dense and drawn from a normal distribution.

    Parameters
    ----------
    n_reservoir : int
        Number of neurons in the reservoir (more = richer features, slower).
    n_genes : int
        Number of input genes (i.e., number of columns in the expression matrix X).
    spectral_radius : float, default=0.9
        Target spectral radius of W_res. Must be < 1 for stable dynamics.
        Values close to 1 give longer memory; values close to 0 forget quickly.
    input_scale : float, default=0.5
        Standard deviation of the W_in entries. Controls how strongly
        gene expression drives the reservoir relative to its own dynamics.
    density : float, default=0.01
        Fraction of non-zero entries in W_res (1% by default = very sparse).
    random_state : int, default=42
        Seed for reproducibility.

    Returns
    -------
    weights : dict with keys:
        'W_res' : scipy.sparse.csr_matrix of shape (n_reservoir, n_reservoir)
            Recurrent reservoir weight matrix.
        'W_in'  : np.ndarray of shape (n_reservoir, n_genes)
            Input weight matrix.
        'n_reservoir' : int
        'n_genes'     : int
        'spectral_radius' : float
        'input_scale'     : float
        'leak_rate'       : float  (stored as 0.3 default; passed to compute_reservoir_states)
        'density'         : float
    """
    rng = np.random.RandomState(random_state)

    # --- Build sparse W_res ---
    W_res = sparse.random(
        n_reservoir, n_reservoir,
        density=density,
        random_state=rng,
        format='csr'
    )

    # Scale so that spectral radius equals the requested value.
    # We estimate ρ(W_res) via the largest singular value of W_res:
    #   σ_max(W_res)² = largest eigenvalue of W_res^T W_res
    eigenvalues = sparse.linalg.eigsh(
        W_res.T @ W_res,
        k=1,
        which='LM',
        return_eigenvectors=False
    )
    rho = np.sqrt(eigenvalues[0])
    W_res = W_res * (spectral_radius / (rho + 1e-10))

    # --- Build dense W_in ---
    # Each gene independently connects to every reservoir neuron
    W_in = rng.randn(n_reservoir, n_genes) * input_scale

    return {
        'W_res':          W_res,
        'W_in':           W_in,
        'n_reservoir':    n_reservoir,
        'n_genes':        n_genes,
        'spectral_radius': spectral_radius,
        'input_scale':    input_scale,
        'density':        density,
    }


# ---------------------------------------------------------------------------
# Core reservoir dynamics  (produces the H matrix)
# ---------------------------------------------------------------------------

def compute_reservoir_states(X, weights, leak_rate=0.3, washout=100):
    """
    Run Echo State Network dynamics over pseudotime-ordered cells.

    At each time step t the reservoir state h(t) is updated by:

        h(t+1) = (1 - α) · h(t)  +  α · tanh( W_res · h(t)  +  W_in · x(t) )

    where α = leak_rate.  This combines:
      - (1 - α) · h(t)          : memory of the previous state
      - α · tanh( W_res · h(t) ): recurrent contribution (echoes of past inputs)
      - α · tanh( W_in · x(t) ) : drive from current cell's gene expression

    The first `washout` steps are discarded so the reservoir "forgets" the
    arbitrary all-zeros initial state before we start recording.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.  Rows must already be sorted by pseudotime
        before calling this function (use order_by_pseudotime from utils).
    weights : dict
        Output of build_reservoir().  Must contain 'W_res' and 'W_in'.
    leak_rate : float, default=0.3
        α in the update equation.  0 = no update (full memory),
        1 = full update (no memory carry-over).
    washout : int, default=100
        Number of initial steps to run before recording states.
        During washout, inputs are cycled (X[t % n_cells]) so the reservoir
        reaches a stable "echo" of the data distribution.

    Returns
    -------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix.
        H[t] = internal state of the reservoir after processing cells 0..t.
        Because of the recurrent dynamics, H[t] encodes a nonlinear summary
        of the entire developmental history up to cell t — not just cell t alone.
    """
    W_res = weights['W_res']
    W_in  = weights['W_in']
    n_reservoir = weights['n_reservoir']

    n_cells = X.shape[0]
    H = np.zeros((n_cells, n_reservoir))   # output: one row per cell
    h = np.zeros(n_reservoir)              # current reservoir state, starts at 0

    total_steps = n_cells + washout

    for t in range(total_steps):

        # Pick the input for this step
        if t < washout:
            x_t = X[t % n_cells]   # cycle through data during washout
        else:
            x_t = X[t - washout]   # actual cell in pseudotime order

        # --- Update equation ---
        drive = W_res @ h + W_in @ x_t          # recurrent + input drive
        h = (1.0 - leak_rate) * h + leak_rate * np.tanh(drive)

        # Record only after washout
        if t >= washout:
            H[t - washout] = h

    return H


# ---------------------------------------------------------------------------
# Graph-regularized reservoir states
# ---------------------------------------------------------------------------

def build_knn_graph(X, n_neighbors=10):
    """
    Build a k-nearest-neighbor graph from a gene expression matrix.

    Cells are connected to their `n_neighbors` nearest neighbors in
    expression space.  Edge weights decay with distance: w = exp(-d²).

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (does NOT need to be pseudotime-sorted).
    n_neighbors : int, default=10
        Number of nearest neighbors per cell.

    Returns
    -------
    adjacency : scipy.sparse.csr_matrix of shape (n_cells, n_cells)
        Weighted adjacency matrix of the k-NN graph.
        adjacency[i, j] = exp(-||X[i] - X[j]||²)  if j is a neighbor of i,
        else 0.
    """
    from sklearn.neighbors import NearestNeighbors

    n_cells = X.shape[0]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)

    row_idx, col_idx, data = [], [], []
    for i in range(n_cells):
        for j, d in zip(indices[i], distances[i]):
            row_idx.append(i)
            col_idx.append(j)
            data.append(np.exp(-d ** 2))

    adjacency = csr_matrix(
        (data, (row_idx, col_idx)), shape=(n_cells, n_cells)
    )
    return adjacency


def compute_graph_reservoir_states(X, weights, adjacency, leak_rate=0.3, washout=100, n_graph_iter=5):
    """
    Compute reservoir states with graph-neighbor smoothing.

    This is a two-step process:
      1. Run standard reservoir dynamics to get initial H.
      2. Iteratively smooth H across the k-NN graph:
             H  ←  tanh( A_norm · H )
         repeated n_graph_iter times.
         A_norm is the row-normalized adjacency matrix (each row sums to 1).

    The effect: each cell's reservoir state is "blurred" toward its
    transcriptionally similar neighbors, providing a form of spatial
    regularization that reduces noise from the pseudotime ordering.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (pseudotime-sorted).
    weights : dict
        Output of build_reservoir().
    adjacency : scipy.sparse.csr_matrix of shape (n_cells, n_cells)
        Output of build_knn_graph().
    leak_rate : float, default=0.3
        Passed to compute_reservoir_states.
    washout : int, default=100
        Passed to compute_reservoir_states.
    n_graph_iter : int, default=5
        Number of graph-smoothing iterations.

    Returns
    -------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Graph-regularized reservoir state matrix.
    """
    # Step 1: standard reservoir pass
    H = compute_reservoir_states(X, weights, leak_rate=leak_rate, washout=washout)

    # Step 2: row-normalize the adjacency matrix
    A = adjacency.astype(np.float32)
    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0                   # avoid division by zero
    A_norm = A.multiply(1.0 / row_sums[:, None])    # row-normalized

    # Step 3: iterative graph smoothing
    for _ in range(n_graph_iter):
        H = np.tanh(A_norm @ H)

    return H
