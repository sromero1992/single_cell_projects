"""
Attractor landscape module for scReservoir.

Models developmental landscapes through identification of attractor states,
energy landscape computation, and fate probability estimation.
"""

from typing import Optional, Tuple
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
from sklearn.utils.extmath import randomized_svd


class ScReservoirLandscape:
    """
    Attractor landscape modeling for single-cell development.

    Estimates latent dynamics, identifies attractor states,
    computes energy landscape, and predicts cell fate probabilities.

    Parameters
    ----------
    n_attractors : int, default=5
        Number of attractor states to identify.
    svd_rank : int, default=50
        Rank for randomized SVD compression of reservoir states.
    lambda_reg : float, default=1e-3
        Regularization parameter for dynamics estimation.
    random_state : int, default=42
        Random seed for reproducibility.
    """

    def __init__(
        self,
        n_attractors: int = 5,
        svd_rank: int = 50,
        lambda_reg: float = 1e-3,
        random_state: int = 42
    ):
        self.n_attractors = n_attractors
        self.svd_rank = svd_rank
        self.lambda_reg = lambda_reg
        self.random_state = random_state

        self.H_latent = None
        self.U = None
        self.S = None
        self.Vt = None
        self.A = None
        self.attractors = None
        self.attractor_cells = None
        self.energy = None
        self.velocity = None
        self.fate_probs = None
        self.pseudotime = None

    def fit(
        self,
        X: np.ndarray,
        H: np.ndarray,
        pseudotime: np.ndarray,
        velocity: Optional[np.ndarray] = None
    ) -> 'ScReservoirLandscape':
        """
        Fit attractor landscape model.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        H : np.ndarray of shape (n_cells, n_reservoir)
            Reservoir state matrix.
        pseudotime : np.ndarray of shape (n_cells,)
            Pseudotime values for cell ordering.
        velocity : np.ndarray of shape (n_cells, n_genes), optional
            Gene velocity. If None, computed from finite differences.

        Returns
        -------
        self
        """
        self.pseudotime = pseudotime.copy()

        # Sort by pseudotime
        sort_idx = np.argsort(pseudotime)
        X_sorted = X[sort_idx]
        H_sorted = H[sort_idx]
        pt_sorted = pseudotime[sort_idx]

        # Compress reservoir states via randomized SVD
        self._compress_states(H_sorted)

        # Compute velocities if not provided
        if velocity is None:
            velocity = self._compute_velocity(X_sorted, pt_sorted)
        else:
            velocity = velocity[sort_idx]

        self.velocity = velocity

        # Estimate latent dynamics
        self._estimate_dynamics(self.H_latent, pt_sorted)

        # Detect attractors
        self._detect_attractors(X_sorted, velocity)

        # Compute energy landscape
        self._compute_energy(self.H_latent)

        # Compute fate probabilities
        self._compute_fate_probabilities()

        return self

    def _compress_states(self, H: np.ndarray) -> None:
        """
        Compress reservoir states using randomized SVD.

        Parameters
        ----------
        H : np.ndarray of shape (n_cells, n_reservoir)
            Reservoir state matrix.
        """
        n_components = min(self.svd_rank, H.shape[0], H.shape[1])

        self.U, self.S, self.Vt = randomized_svd(
            H, n_components=n_components, random_state=self.random_state
        )

        # Project to latent space
        self.H_latent = self.U @ np.diag(self.S)

    def _compute_velocity(
        self,
        X: np.ndarray,
        pseudotime: np.ndarray
    ) -> np.ndarray:
        """
        Compute gene velocity via finite differences.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Sorted gene expression matrix.
        pseudotime : np.ndarray of shape (n_cells,)
            Sorted pseudotime.

        Returns
        -------
        velocity : np.ndarray of shape (n_cells, n_genes)
            Gene velocity (dX/dt).
        """
        dt = np.diff(pseudotime)
        dX = np.diff(X, axis=0)
        velocity = np.zeros_like(X)

        velocity[1:] = dX / dt[:, None]
        velocity[0] = velocity[1]

        return velocity

    def _estimate_dynamics(
        self,
        H_latent: np.ndarray,
        pseudotime: np.ndarray
    ) -> None:
        """
        Estimate latent dynamics matrix A from compressed states.

        Parameters
        ----------
        H_latent : np.ndarray of shape (n_cells, svd_rank)
            Compressed reservoir states.
        pseudotime : np.ndarray of shape (n_cells,)
            Pseudotime values.
        """
        # dH/dt = H @ A + noise
        # Using finite differences: (H_{t+1} - H_t) / dt = H_t @ A

        dt = np.diff(pseudotime)
        dH = np.diff(H_latent, axis=0)

        # Weighted least squares: minimize ||dH - H[:-1] @ A||^2 + lambda*||A||^2
        H_prev = H_latent[:-1]

        # Solve: (H.T @ H + lambda*I) @ A.T = H.T @ dH
        HTH = H_prev.T @ H_prev
        HTdH = H_prev.T @ dH / dt[:, None]

        A_T = np.linalg.solve(
            HTH + self.lambda_reg * np.eye(HTH.shape[0]),
            HTdH
        )

        self.A = A_T.T

    def _detect_attractors(
        self,
        X: np.ndarray,
        velocity: np.ndarray
    ) -> None:
        """
        Detect attractor states as low-velocity cells.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Gene expression matrix.
        velocity : np.ndarray of shape (n_cells, n_genes)
            Gene velocities.
        """
        # Compute velocity magnitude
        vel_mag = np.linalg.norm(velocity, axis=1)

        # Find cells with low velocity (potential attractors)
        vel_threshold = np.percentile(vel_mag, 10)  # Bottom 10%
        attractor_mask = vel_mag < vel_threshold

        attractor_cells = np.where(attractor_mask)[0]

        # Cluster attractor cells to find attractor states
        if len(attractor_cells) > self.n_attractors:
            X_attractors = X[attractor_cells]
            kmeans = KMeans(
                n_clusters=self.n_attractors,
                random_state=self.random_state,
                n_init=10
            )
            clusters = kmeans.fit_predict(X_attractors)

            # Store attractor states as cluster centers
            self.attractors = kmeans.cluster_centers_
            self.attractor_cells = [
                attractor_cells[clusters == i]
                for i in range(self.n_attractors)
            ]
        else:
            # Too few attractor cells
            self.attractors = X[attractor_cells]
            self.attractor_cells = [attractor_cells]

    def _compute_energy(self, H_latent: np.ndarray) -> None:
        """
        Compute energy landscape via kernel density estimation.

        Parameters
        ----------
        H_latent : np.ndarray of shape (n_cells, svd_rank)
            Compressed reservoir states.
        """
        if H_latent.shape[0] < 10:
            self.energy = np.zeros(H_latent.shape[0])
            return

        # Use KDE to estimate probability density
        try:
            kde = gaussian_kde(H_latent.T)
            energy = -np.log(kde(H_latent.T) + 1e-10)
        except np.linalg.LinAlgError:
            # If KDE fails, use Euclidean distance to attractor centers
            energy = np.zeros(H_latent.shape[0])
            for i, attractor in enumerate(self.attractors):
                dist = np.linalg.norm(H_latent - attractor[np.newaxis, :], axis=1)
                energy += dist

        self.energy = energy

    def _compute_fate_probabilities(self) -> None:
        """
        Compute cell fate probabilities based on distance to attractors.
        """
        n_cells = self.H_latent.shape[0]
        self.fate_probs = np.zeros((n_cells, len(self.attractors)))

        for i, attractor in enumerate(self.attractors):
            # Euclidean distance to attractor
            dist = np.linalg.norm(
                self.H_latent - attractor[np.newaxis, :], axis=1
            )

            # Convert to probability via softmax
            self.fate_probs[:, i] = np.exp(-dist) / (1 + np.exp(-dist))

        # Normalize across attractors
        row_sum = self.fate_probs.sum(axis=1, keepdims=True)
        self.fate_probs = self.fate_probs / (row_sum + 1e-10)

    def get_energy_landscape(self) -> np.ndarray:
        """
        Get energy landscape (low values = attractors).

        Returns
        -------
        energy : np.ndarray of shape (n_cells,)
            Energy values per cell.
        """
        if self.energy is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.energy.copy()

    def get_fate_probabilities(self) -> np.ndarray:
        """
        Get cell fate probabilities.

        Returns
        -------
        fate_probs : np.ndarray of shape (n_cells, n_attractors)
            Probability of each cell committing to each attractor.
        """
        if self.fate_probs is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.fate_probs.copy()

    def get_latent_states(self) -> np.ndarray:
        """
        Get compressed latent states.

        Returns
        -------
        H_latent : np.ndarray of shape (n_cells, svd_rank)
            Compressed reservoir states.
        """
        if self.H_latent is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.H_latent.copy()

    def get_dynamics_matrix(self) -> np.ndarray:
        """
        Get estimated latent dynamics matrix.

        Returns
        -------
        A : np.ndarray of shape (svd_rank, svd_rank)
            Dynamics matrix where dH/dt ≈ H @ A.
        """
        if self.A is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.A.copy()

    def get_attractor_states(self) -> np.ndarray:
        """
        Get attractor gene expression states.

        Returns
        -------
        attractors : np.ndarray of shape (n_attractors, n_genes)
            Gene expression profiles of attractors.
        """
        if self.attractors is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")
        return self.attractors.copy()

    def get_attractor_gene_programs(
        self,
        X: np.ndarray
    ) -> np.ndarray:
        """
        Get mean gene expression of cells in each attractor basin.

        Parameters
        ----------
        X : np.ndarray of shape (n_cells, n_genes)
            Original gene expression matrix.

        Returns
        -------
        programs : np.ndarray of shape (n_attractors, n_genes)
            Mean gene expression per attractor.
        """
        if self.attractor_cells is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        n_genes = X.shape[1]
        programs = np.zeros((len(self.attractors), n_genes))

        for i, cells in enumerate(self.attractor_cells):
            if len(cells) > 0:
                programs[i] = X[cells].mean(axis=0)

        return programs

    def get_grn_from_dynamics(self, W_in: np.ndarray) -> np.ndarray:
        """
        Infer GRN from dynamics matrix.

        GRN = |W_in.T @ A @ W_in| represents gene regulatory relationships
        encoded in the latent dynamics.

        Parameters
        ----------
        W_in : np.ndarray of shape (n_reservoir, n_genes)
            Input weight matrix from reservoir.

        Returns
        -------
        GRN : np.ndarray of shape (n_genes, n_genes)
            Gene influence matrix.
        """
        if self.A is None:
            raise RuntimeError("Model not fitted yet. Call fit() first.")

        # Project dynamics through input/output weights
        # Gene i -> Reservoir -> Latent -> Gene j
        if self.U is None:
            raise RuntimeError("Latent compression not performed.")

        # Full dynamics in reservoir space: W_in.T @ A @ W_in
        # But we have compressed dynamics, so use U to project back
        A_expanded = self.U @ self.A @ self.U.T

        GRN = np.abs(W_in.T @ A_expanded @ W_in)

        # Normalize
        GRN_max = np.max(GRN)
        if GRN_max > 0:
            GRN = GRN / GRN_max

        return GRN
