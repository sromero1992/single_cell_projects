"""
Gene Regulatory Network (GRN) inference functions for scReservoir.

All functions are stateless: inputs in, outputs out.
The central idea is a two-step process:
  1. fit_readout_*()  — learn W_out by regressing gene expression on H
  2. infer_grn()      — back-project W_out through W_in to get gene-gene influences

Three regression modes are available:
  - 'standard' : regress X directly on H  (no temporal structure assumed)
  - 'causal'   : regress X[t] on H[t-1]  (causal / Granger-style)
  - 'velocity' : regress dX/dt on H       (RNA-velocity style)
"""

import numpy as np
from utils import ridge_regression, normalize_grn, order_by_pseudotime


# ---------------------------------------------------------------------------
# Readout weight fitting  (the only "training" in reservoir computing)
# ---------------------------------------------------------------------------

def fit_readout_standard(H, X, lambda_reg=1e-3):
    """
    Fit output weights by regressing gene expression directly on reservoir states.

    Solves:  minimize  ||H @ W_out - X||²  +  lambda * ||W_out||²

    No temporal structure is used here: each cell is treated independently.
    Use this when you do not have pseudotime or do not need causal inference.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix from compute_reservoir_states().
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (same ordering as H).
    lambda_reg : float, default=1e-3
        Ridge regularization strength.

    Returns
    -------
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Learned readout weight matrix.
    """
    return ridge_regression(H, X, lambda_reg)


def fit_readout_causal(H, X, pseudotime, lambda_reg=1e-3):
    """
    Fit output weights using time-lagged (causal) regression.

    Regresses the NEXT cell's expression on the CURRENT reservoir state:
        X[t]  ~  H[t-1] @ W_out

    By using the previous reservoir state to predict the current expression,
    we enforce a temporal arrow of causality: gene A can only influence gene B
    if A's signal in the reservoir (H) precedes B's expression.
    This is a Granger-causality approach adapted for pseudotime.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir states.  Will be sorted by pseudotime internally.
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime per cell.  Used to establish the temporal ordering.
    lambda_reg : float, default=1e-3
        Ridge regularization strength.

    Returns
    -------
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Learned readout weight matrix.
    """
    # Sort both H and X by pseudotime
    sort_idx = np.argsort(pseudotime)
    H_sorted = H[sort_idx]
    X_sorted = X[sort_idx]

    # Lag: use reservoir state at t-1 to predict expression at t
    H_lagged  = H_sorted[:-1]   # all cells except the last  → "past" states
    X_target  = X_sorted[1:]    # all cells except the first → "future" expression

    return ridge_regression(H_lagged, X_target, lambda_reg)


def fit_readout_velocity(H, X, pseudotime, lambda_reg=1e-3):
    """
    Fit output weights by regressing gene velocity on midpoint reservoir states.

    Regresses dX/dt on the average of adjacent reservoir states:
        dX/dt  ~  H_midpoint @ W_out
    where H_midpoint[t] = (H[t] + H[t+1]) / 2.

    This is the most RNA-velocity-aligned mode: instead of predicting
    expression levels, W_out predicts the *rate of change* of expression.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir states.
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime per cell.
    lambda_reg : float, default=1e-3
        Ridge regularization strength.

    Returns
    -------
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Learned readout weight matrix (predicts velocity, not expression).
    """
    sort_idx  = np.argsort(pseudotime)
    X_sorted  = X[sort_idx]
    H_sorted  = H[sort_idx]
    pt_sorted = pseudotime[sort_idx]

    # Finite-difference velocity: dX/dt between adjacent pseudotime steps
    dt       = np.diff(pt_sorted)                   # shape (n_cells-1,)
    dX       = np.diff(X_sorted, axis=0)            # shape (n_cells-1, n_genes)
    velocity = dX / dt[:, None]                     # element-wise division

    # Midpoint reservoir state between each adjacent pair
    H_mid = (H_sorted[:-1] + H_sorted[1:]) / 2.0   # shape (n_cells-1, n_reservoir)

    return ridge_regression(H_mid, velocity, lambda_reg)


# ---------------------------------------------------------------------------
# GRN back-projection
# ---------------------------------------------------------------------------

def infer_grn(W_in, W_out, normalize=True):
    """
    Infer a gene-gene influence matrix by back-projecting readout weights.

    The logic:
      - W_in[k, i]  = how strongly gene i activates reservoir neuron k
      - W_out[k, j] = how strongly reservoir neuron k predicts gene j
      - Therefore, the influence of gene i on gene j through all reservoir
        neurons is the dot product:  sum_k  W_in[k,i] * W_out[k,j]

    In matrix form:
        GRN  =  |W_in.T  @  W_out|         (n_genes × n_genes)

    GRN[i, j] ≈ influence of gene i on gene j.

    Parameters
    ----------
    W_in : np.ndarray of shape (n_reservoir, n_genes)
        Fixed input weight matrix from build_reservoir().
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Fitted readout weights from any of the fit_readout_* functions.
    normalize : bool, default=True
        If True, apply min-max normalization to scale GRN to [0, 1].

    Returns
    -------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene influence matrix.  GRN[i, j] = influence of gene i on gene j.
    """
    GRN = np.abs(W_in.T @ W_out)   # shape (n_genes, n_genes)

    if normalize:
        GRN = normalize_grn(GRN, method='minmax')

    return GRN


def infer_grn_from_dynamics(W_in, A, U, normalize=True):
    """
    Derive a GRN from the latent dynamics matrix estimated in the landscape module.

    The latent dynamics matrix A encodes how the compressed reservoir state
    evolves over pseudotime.  We project it back to gene space through the
    input weights to get a dynamics-based GRN.

    Math:
        A_expanded  =  U  @  A  @  U.T        (back-project from SVD latent space)
        GRN         =  |W_in.T  @  A_expanded  @  W_in|

    Parameters
    ----------
    W_in : np.ndarray of shape (n_reservoir, n_genes)
        Fixed input weight matrix.
    A : np.ndarray of shape (svd_rank, svd_rank)
        Latent dynamics matrix from estimate_latent_dynamics().
    U : np.ndarray of shape (n_cells, svd_rank)
        Left singular vectors from compress_reservoir_states().
    normalize : bool, default=True
        If True, apply min-max normalization.

    Returns
    -------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Dynamics-based gene influence matrix.
    """
    A_expanded = U @ A @ U.T                        # (n_cells, n_cells)
    GRN = np.abs(W_in.T @ A_expanded @ W_in)        # (n_genes, n_genes)

    if normalize:
        GRN = normalize_grn(GRN, method='minmax')

    return GRN


# ---------------------------------------------------------------------------
# Prediction and evaluation
# ---------------------------------------------------------------------------

def predict_expression(H, W_out):
    """
    Predict gene expression from reservoir states using fitted readout weights.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir states (can be from new, unseen cells).
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Fitted readout weights from any fit_readout_* function.

    Returns
    -------
    X_pred : np.ndarray of shape (n_cells, n_genes)
        Predicted gene expression.
    """
    return H @ W_out


def compute_prediction_mse(X_true, H, W_out):
    """
    Compute mean squared error of gene expression predictions.

    Parameters
    ----------
    X_true : np.ndarray of shape (n_cells, n_genes)
        True gene expression.
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir states corresponding to X_true.
    W_out : np.ndarray of shape (n_reservoir, n_genes)
        Fitted readout weights.

    Returns
    -------
    mse : float
        Mean squared error averaged over all cells and genes.
    """
    X_pred = predict_expression(H, W_out)
    return float(np.mean((X_true - X_pred) ** 2))
