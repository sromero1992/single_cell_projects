"""
Core reservoir computing module for scReservoir.

Implements Echo State Networks with random recurrent connections
and trainable linear readout for single-cell RNA-seq applications.
"""

from typing import Optional, Tuple, Union
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
from sklearn.utils.extmath import randomized_svd


class ScReservoir:
    """
    Reservoir computing module for single-cell RNA-seq data.

    Implements an Echo State Network with random recurrent connections.
    The reservoir computes high-dimensional nonlinear transformations
    of input gene expression data, which are then used for regression.

    Parameters
    ----------
    n_reservoir : int, default=500
        Number of neurons in the recurrent reservoir.
    spectral_radius : float, default=0.9
        Spectral radius of the recurrent weight matrix W_res.
        Controls the stability of dynamics (should be < 1).
    input_scale : float, default=0.5
        Scaling factor for input weight matrix W_in.
    leak_rate : float, default=0.3
        Leak rate for reservoir neurons (0 = no feedback, 1 = full memory).
    density : float, default=0.01
        Sparsity of the recurrent connection matrix (fraction of non-zero weights).
    random_state : int, default=42
        Random seed for reproducibility.
    sparse_mode : bool, default=True
        If True, use sparse matrices for reservoir weights.
    """

    def __init__(
        self,
        n_reservoir: int = 500,
        spectral_radius: float = 0.9,
        input_scale: float = 0.5,
        leak_rate: float = 0.3,
        density: float = 0.01,
        random_state: int = 42,
        sparse_mode: bool = True
    ):
        self.n_reservoir = n_reservoir
        self.spectral_radius = spectral_radius
        self.input_scale = input_scale
        self.leak_rate = leak_rate
        self.density = density
        self.random_state = random_state
        self.sparse_mode = sparse_mode

        self.W_res = None
        self.W_in = None
        self.n_genes = None
        self._rng = np.random.RandomState(random_state)

    def _initialize_weights(self, n_genes: int) -> None:
        """
        Initialize reservoir and input weight matrices.

        Parameters
        ----------
        n_genes : int
            Number of input genes (dimensions).
        """
        self.n_genes = n_genes

        if self.sparse_mode:
            # Generate sparse random matrix
            self.W_res = sparse.random(
                self.n_reservoir, self.n_reservoir,
                density=self.density,
                random_state=self._rng,
                format='csr'
            )
            # Scale to desired spectral radius
            eigenvalues = sparse.linalg.eigsh(
                self.W_res.T @ self.W_res,
                k=1, which='LM', return_eigenvectors=False
            )
            rho = np.sqrt(eigenvalues[0])
            self.W_res *= self.spectral_radius / (rho + 1e-10)
        else:
            # Dense reservoir
            W_res_dense = self._rng.randn(self.n_reservoir, self.n_reservoir)
            # Make sparse
            mask = self._rng.rand(self.n_reservoir, self.n_reservoir) < self.density
            W_res_dense[~mask] = 0
            # Normalize by spectral radius
            eigenvalues = np.linalg.eigvals(W_res_dense)
            rho = np.max(np.abs(eigenvalues))
            W_res_dense = W_res_dense * self.spectral_radius / (rho + 1e-10)
            self.W_res = csr_matrix(W_res_dense) if self.sparse_mode else W_res_dense

        # Input weights
        self.W_in = self._rng.randn(self.n_reservoir, n_genes) * self.input_scale

    def compute_states(
        self,
        X: np.ndarray,
        pseudotime: Optional[np.ndarray] = None,
        washout: int = 100
    ) -> np.ndarray:
        """
        Run reservoir dynamics and compute state matrix.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        pseudotime : np.ndarray of shape (n_cells,), optional
            Pseudotime values for ordering cells. If None, uses cell order.
        washout : int, default=100
            Number of initial time steps to discard (transient behavior).

        Returns
        -------
        H : np.ndarray of shape (n_cells, n_reservoir)
            Reservoir state matrix.
        """
        if self.W_res is None:
            self._initialize_weights(X.shape[1])

        n_cells = X.shape[0]
        H = np.zeros((n_cells, self.n_reservoir))
        h = np.zeros(self.n_reservoir)

        for t in range(n_cells + washout):
            # Current input
            if t < washout:
                # Washout phase: cycle through data
                x_t = X[t % n_cells]
            else:
                x_t = X[t - washout]

            # Reservoir update: h(t+1) = (1-leak)*h(t) + leak*tanh(W_res@h(t) + W_in@x(t))
            if isinstance(self.W_res, csr_matrix):
                h_new = self.W_res @ h
            else:
                h_new = self.W_res @ h

            h_new = h_new + self.W_in @ x_t
            h = (1 - self.leak_rate) * h + self.leak_rate * np.tanh(h_new)

            if t >= washout:
                H[t - washout] = h

        return H

    def compute_states_graph(
        self,
        X: np.ndarray,
        adjacency: Union[np.ndarray, csr_matrix],
        n_iter: int = 5
    ) -> np.ndarray:
        """
        Compute reservoir states using graph propagation.

        Uses k-nearest neighbor graph to propagate information across cells.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        adjacency : np.ndarray or csr_matrix of shape (n_cells, n_cells)
            Adjacency matrix of cell neighborhood graph (typically k-NN).
        n_iter : int, default=5
            Number of graph propagation iterations.

        Returns
        -------
        H : np.ndarray of shape (n_cells, n_reservoir)
            Graph-regularized reservoir state matrix.
        """
        if self.W_res is None:
            self._initialize_weights(X.shape[1])

        n_cells = X.shape[0]
        H = np.zeros((n_cells, self.n_reservoir))

        # Initial reservoir computation
        h_init = self.compute_states(X)
        H = h_init.copy()

        # Graph propagation
        A_norm = adjacency.astype(np.float32)
        if isinstance(A_norm, csr_matrix):
            row_sum = np.array(A_norm.sum(axis=1)).flatten()
        else:
            row_sum = A_norm.sum(axis=1)

        row_sum[row_sum == 0] = 1
        if isinstance(A_norm, csr_matrix):
            A_norm = A_norm.multiply(1.0 / row_sum[:, None])
        else:
            A_norm = A_norm / row_sum[:, None]

        for _ in range(n_iter):
            if isinstance(A_norm, csr_matrix):
                H = A_norm @ H
            else:
                H = A_norm @ H
            H = np.tanh(H)  # Nonlinearity

        return H


class ScGraphReservoir(ScReservoir):
    """
    Graph-enhanced reservoir computing for single-cell data.

    Extends ScReservoir to incorporate cell neighborhood information
    through graph-based propagation of reservoir states.

    Inherits all parameters from ScReservoir.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.adjacency = None

    def fit_graph(
        self,
        X: np.ndarray,
        n_neighbors: int = 10
    ) -> None:
        """
        Build k-nearest neighbor graph from expression data.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        n_neighbors : int, default=10
            Number of nearest neighbors for each cell.
        """
        from sklearn.neighbors import NearestNeighbors

        nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(X)
        distances, indices = nbrs.kneighbors(X)

        # Build adjacency matrix
        n_cells = X.shape[0]
        row_idx = []
        col_idx = []
        data = []

        for i in range(n_cells):
            for j, d in zip(indices[i], distances[i]):
                weight = np.exp(-d**2)
                row_idx.append(i)
                col_idx.append(j)
                data.append(weight)

        self.adjacency = csr_matrix(
            (data, (row_idx, col_idx)), shape=(n_cells, n_cells)
        )

    def compute_states(
        self,
        X: np.ndarray,
        pseudotime: Optional[np.ndarray] = None,
        use_graph: bool = True,
        n_graph_iter: int = 5
    ) -> np.ndarray:
        """
        Compute reservoir states with optional graph propagation.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        pseudotime : np.ndarray of shape (n_cells,), optional
            Pseudotime values.
        use_graph : bool, default=True
            If True and adjacency exists, use graph propagation.
        n_graph_iter : int, default=5
            Number of graph propagation iterations.

        Returns
        -------
        H : np.ndarray of shape (n_cells, n_reservoir)
            Reservoir state matrix.
        """
        if use_graph and self.adjacency is not None:
            return self.compute_states_graph(X, self.adjacency, n_iter=n_graph_iter)
        else:
            return super().compute_states(X, pseudotime)
