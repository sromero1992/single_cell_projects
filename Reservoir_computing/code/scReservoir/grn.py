"""
Gene Regulatory Network (GRN) inference module for scReservoir.

Implements GRN inference through regression of gene expression on
reservoir states, with support for standard, causal, and velocity modes.
"""

from typing import Optional, Dict, Tuple
import numpy as np
from .utils import ridge_regression, normalize_grn


class ScReservoirGRN:
    """
    Gene Regulatory Network inference using reservoir computing.

    Infers GRN by regressing gene expression on reservoir states
    and back-projecting readout weights through input connections.

    Parameters
    ----------
    reservoir : ScReservoir
        Fitted reservoir computing model.
    lambda_reg : float, default=1e-3
        L2 regularization parameter for ridge regression.
    n_top_regulators : int, default=10
        Number of top regulators to return in summary statistics.
    """

    def __init__(
        self,
        reservoir,
        lambda_reg: float = 1e-3,
        n_top_regulators: int = 10
    ):
        self.reservoir = reservoir
        self.lambda_reg = lambda_reg
        self.n_top_regulators = n_top_regulators

        self.W_out = None  # Readout weights (n_reservoir, n_genes)
        self.GRN = None    # Gene influence matrix (n_genes, n_genes)
        self.H = None      # Reservoir states
        self.X = None      # Expression data
        self.mode = None   # Fitting mode used

    def fit(
        self,
        X: np.ndarray,
        H: np.ndarray,
        pseudotime: Optional[np.ndarray] = None,
        mode: str = 'standard'
    ) -> 'ScReservoirGRN':
        """
        Fit GRN by regressing expression on reservoir states.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        H : np.ndarray of shape (n_cells, n_reservoir)
            Reservoir state matrix (from reservoir.compute_states).
        pseudotime : np.ndarray of shape (n_cells,), optional
            Pseudotime values. Required for 'causal' mode.
        mode : {'standard', 'causal', 'velocity'}, default='standard'
            Fitting mode:
            - 'standard': regress X on H
            - 'causal': regress X on time-lagged H (enforce causality)
            - 'velocity': regress dX/dt on H at midpoints

        Returns
        -------
        self
        """
        self.X = X
        self.H = H
        self.mode = mode
        n_cells, n_genes = X.shape

        if mode == 'standard':
            # Standard regression: X = H @ W_out + noise
            self.W_out = ridge_regression(H, X, self.lambda_reg)

        elif mode == 'causal':
            if pseudotime is None:
                raise ValueError("pseudotime required for causal mode")

            # Sort by pseudotime
            sort_idx = np.argsort(pseudotime)
            X_sorted = X[sort_idx]
            H_sorted = H[sort_idx]

            # Lag reservoir states: H_t-1 -> X_t
            # Skip first cell (no previous state)
            H_lagged = H_sorted[:-1]
            X_target = X_sorted[1:]

            self.W_out = ridge_regression(H_lagged, X_target, self.lambda_reg)

        elif mode == 'velocity':
            # Velocity-based: dX/dt = H @ W_out
            if pseudotime is None:
                raise ValueError("pseudotime required for velocity mode")

            # Compute finite difference velocities
            sort_idx = np.argsort(pseudotime)
            X_sorted = X[sort_idx]
            pt_sorted = pseudotime[sort_idx]
            H_sorted = H[sort_idx]

            # Compute dX/dt
            dt = np.diff(pt_sorted)
            dX = np.diff(X_sorted, axis=0)
            velocity = dX / dt[:, None]

            # Use midpoint reservoir states
            H_mid = (H_sorted[:-1] + H_sorted[1:]) / 2

            self.W_out = ridge_regression(H_mid, velocity, self.lambda_reg)

        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Back-project to get gene-gene influence
        self._compute_grn()

        return self

    def _compute_grn(self) -> None:
        """
        Back-project readout weights through input connections to infer GRN.

        GRN_ij approximates influence of gene i on gene j through
        the reservoir: GRN = |W_in.T @ W_out|
        """
        if self.W_out is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        # Back-projection: GRN = |W_in.T @ W_out|
        W_in = self.reservoir.W_in  # shape (n_reservoir, n_genes)
        W_out = self.W_out  # shape (n_reservoir, n_genes)

        # Compute influence: gene_i -> reservoir -> gene_j
        # GRN[i,j] = |sum_k W_in[k,i] * W_out[k,j]|
        GRN = np.abs(W_in.T @ W_out)  # shape (n_genes, n_genes)

        # Normalize
        self.GRN = normalize_grn(GRN)

    def get_grn(self) -> np.ndarray:
        """
        Get inferred gene regulatory network matrix.

        Returns
        -------
        GRN : np.ndarray of shape (n_genes, n_genes)
            Gene regulatory network matrix, normalized to [0, 1].
            GRN[i, j] represents the influence of gene i on gene j.
        """
        if self.GRN is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.GRN.copy()

    def get_top_regulators(
        self,
        gene_names: Optional[np.ndarray] = None,
        k: Optional[int] = None
    ) -> Dict[str, np.ndarray]:
        """
        Get top k regulators for each gene.

        Parameters
        ----------
        gene_names : np.ndarray of shape (n_genes,), optional
            Gene names/IDs. If None, uses integer indices.
        k : int, optional
            Number of top regulators. Defaults to n_top_regulators.

        Returns
        -------
        regulators : dict
            Dictionary mapping gene names to array of (regulator, influence_score) tuples.
        """
        if self.GRN is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        k = k or self.n_top_regulators
        n_genes = self.GRN.shape[0]

        if gene_names is None:
            gene_names = np.arange(n_genes).astype(str)

        regulators = {}

        for j in range(n_genes):
            target_gene = gene_names[j]
            influences = self.GRN[:, j]
            top_idx = np.argsort(influences)[-k:][::-1]

            regulators[target_gene] = [
                (gene_names[i], influences[i]) for i in top_idx
            ]

        return regulators

    def get_grn_sparse(self, threshold: float = 0.1) -> np.ndarray:
        """
        Get thresholded (sparse) GRN.

        Parameters
        ----------
        threshold : float, default=0.1
            Influence threshold. Values below threshold are set to 0.

        Returns
        -------
        GRN_sparse : np.ndarray of shape (n_genes, n_genes)
            Sparse GRN with weak interactions removed.
        """
        if self.GRN is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        GRN_sparse = self.GRN.copy()
        GRN_sparse[GRN_sparse < threshold] = 0

        return GRN_sparse

    def get_readout_weights(self) -> np.ndarray:
        """
        Get raw readout weights from reservoir to genes.

        Returns
        -------
        W_out : np.ndarray of shape (n_reservoir, n_genes)
            Linear readout weights.
        """
        if self.W_out is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.W_out.copy()

    def get_input_weights(self) -> np.ndarray:
        """
        Get input weight matrix.

        Returns
        -------
        W_in : np.ndarray of shape (n_reservoir, n_genes)
            Input weights from genes to reservoir.
        """
        return self.reservoir.W_in.copy()

    def predict(self, H_new: np.ndarray) -> np.ndarray:
        """
        Predict gene expression from new reservoir states.

        Parameters
        ----------
        H_new : np.ndarray of shape (n_cells_new, n_reservoir)
            Reservoir states for new cells.

        Returns
        -------
        X_pred : np.ndarray of shape (n_cells_new, n_genes)
            Predicted gene expression.
        """
        if self.W_out is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        return H_new @ self.W_out

    def compute_prediction_error(
        self,
        X_test: np.ndarray,
        H_test: np.ndarray
    ) -> float:
        """
        Compute mean squared error of predictions on test data.

        Parameters
        ----------
        X_test : np.ndarray of shape (n_test, n_genes)
            Test gene expression.
        H_test : np.ndarray of shape (n_test, n_reservoir)
            Test reservoir states.

        Returns
        -------
        mse : float
            Mean squared prediction error.
        """
        X_pred = self.predict(H_test)
        mse = np.mean((X_test - X_pred) ** 2)
        return mse
