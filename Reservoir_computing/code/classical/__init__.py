"""
scReservoir — Classical Reservoir Computing (Echo State Network) module.

All functions are plain Python functions — no classes.

Typical usage
-------------
import sys; sys.path.insert(0, 'path/to/code/classical')
import utils, reservoir, grn, landscape, plotting

# 1. Preprocess
X = utils.preprocess(adata)

# 2. Sort by pseudotime
X_s, pt_s, sort_idx = utils.order_by_pseudotime(X, pseudotime)

# 3. Build reservoir weights (fixed, never trained)
weights = reservoir.build_reservoir(n_reservoir=500, n_genes=X.shape[1])

# 4. Compute H matrix (reservoir state time series)
H = reservoir.compute_reservoir_states(X_s, weights)

# 5. Fit readout (causal mode for developmental data)
W_out = grn.fit_readout_causal(H, X_s, pt_s)

# 6. Infer GRN
GRN = grn.infer_grn(weights['W_in'], W_out)

# 7. Get top regulators
top_regs = utils.get_top_regulators(GRN, gene_names, k=10)

# 8. Full landscape analysis
result = landscape.run_landscape_pipeline(X_s, H, pt_s)
"""
