import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.stats import norm

def sc_norm(X, norm_type='libsize'):
    """Normalize count matrix by library size."""
    if norm_type == 'libsize':
        lib_size = np.sum(X, axis=0)
        lib_size[lib_size == 0] = 1
        X_norm = X / lib_size[np.newaxis, :] * np.median(lib_size)
        return X_norm
    else:
        raise ValueError(f"Normalization type '{norm_type}' not implemented")

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    ranks = np.arange(1, n + 1)
    adjusted = sorted_pvals * n / ranks
    adjusted = np.minimum(adjusted, 1.0)
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i+1])
    results = np.zeros(n)
    results[sorted_indices] = adjusted
    return results

def sc_splinefit(X, genes, sortit=True, plotit=False):
    """
    Fit cubic spline and calculate Per-Gene statistics (d, pval, fdr).
    """
    mean_expr = np.mean(X, axis=1)
    # Add epsilon to avoid division by zero
    cv = np.std(X, axis=1) / (mean_expr + 1e-10)
    dropout_rate = np.sum(X == 0, axis=1) / X.shape[1]
    
    lgu = np.log1p(mean_expr)
    lgcv = np.log1p(cv)
    
    T = pd.DataFrame({
        'genes': genes,
        'lgu': lgu,
        'lgcv': lgcv,
        'dropr': dropout_rate,
        'mean_expr': mean_expr,
        'cv': cv
    })
    
    if sortit:
        T = T.sort_values('lgu').reset_index(drop=True)
        X = X[T.index.values, :]
        genes = T['genes'].values
    
    # 3D Coordinates for spline fitting
    xyz = np.column_stack([T['lgu'], T['lgcv'], T['dropr']])
    t = np.arange(len(T))
    
    # Fit Splines (smoothing factor s is critical)
    # Set A likely used s = len * 0.01 based on typical defaults
    s_factor = len(t) * 0.01
    spline_x = interpolate.UnivariateSpline(t, xyz[:, 0], k=3, s=s_factor)
    spline_y = interpolate.UnivariateSpline(t, xyz[:, 1], k=3, s=s_factor)
    spline_z = interpolate.UnivariateSpline(t, xyz[:, 2], k=3, s=s_factor)
    
    t_fit = np.linspace(0, len(T)-1, len(T))
    xyzp = np.column_stack([spline_x(t_fit), spline_y(t_fit), spline_z(t_fit)])
    
    # Calculate Distances (d)
    distances = np.zeros(len(T))
    nearidx = np.zeros(len(T), dtype=int)
    
    for i in range(len(T)):
        dists = np.linalg.norm(xyzp - xyz[i], axis=1)
        nearidx[i] = np.argmin(dists)
        distances[i] = dists[nearidx[i]]
    
    # --- STATISTICAL FIX: Calculate P-values for deviation from spline ---
    # Z-score of the distance
    dist_z = (distances - np.mean(distances)) / np.std(distances)
    # One-tailed p-value (testing if gene is "noisier" than expected)
    pvals = 1 - norm.cdf(dist_z)
    fdrs = bh_adjust(pvals)
    
    T['d'] = distances       # matches 'd_CRC' in Set A
    T['pval'] = pvals        # matches 'pval_CRC' in Set A
    T['fdr'] = fdrs          # matches 'fdr_CRC' in Set A
    T['nearidx'] = nearidx
    
    return T, X, genes, xyzp

def sc_splinedv(X1, X2, g1, g2, name1='Set1', name2='Set2', fname=None, plotit=False):
    """
    Compute differential variability with full Set A statistics.
    """
    g1 = np.array(g1)
    g2 = np.array(g2)
    g3 = np.intersect1d(g1, g2)
    
    irows = np.where(np.isin(g1, g3))[0]
    jrows = np.where(np.isin(g2, g3))[0]
    
    X1 = X1[irows, :]
    X2 = X2[jrows, :]
    
    X1 = sc_norm(X1, 'libsize')
    X2 = sc_norm(X2, 'libsize')
    
    # Fit Set 1
    T1, X1, _, xyzp1 = sc_splinefit(X1, g3)
    T1 = T1.sort_values('genes').reset_index(drop=True)
    xyz1 = T1[['lgu', 'lgcv', 'dropr']].values
    
    # Fit Set 2
    T2, X2, _, xyzp2 = sc_splinefit(X2, g3)
    T2 = T2.sort_values('genes').reset_index(drop=True)
    xyz2 = T2[['lgu', 'lgcv', 'dropr']].values
    
    # Differential Calculation
    v1 = xyz1 - xyzp1[T1['nearidx'].values, :]
    v2 = xyz2 - xyzp2[T2['nearidx'].values, :]
    
    dist_diff = np.linalg.norm(v1 - v2, axis=1)
    dist_sign = np.sign(np.linalg.norm(v1, axis=1) - np.linalg.norm(v2, axis=1))
    
    # P-value for the DIFFERENCE
    ddz = (dist_diff - np.mean(dist_diff)) / np.std(dist_diff)
    pval_diff = 1 - norm.cdf(ddz)
    
    # Rename Columns
    T1_renamed = T1.add_suffix(f'_{name1}')
    T2_renamed = T2.add_suffix(f'_{name2}')
    
    # Combine
    Tdv = pd.concat([
        T1_renamed,
        T2_renamed,
        pd.DataFrame({'DiffDist': dist_diff, 'DiffSign': dist_sign, 'pval': pval_diff})
    ], axis=1)
    
    # Zero out endpoints (artifacts)
    max_idx1 = T1['nearidx'].max()
    max_idx2 = T2['nearidx'].max()
    idxx = (T1['nearidx'] < 5) | (T2['nearidx'] < 5) | \
           (T1['nearidx'] > max_idx1-5) | (T2['nearidx'] > max_idx2-5)
    Tdv.loc[idxx, 'DiffDist'] = 0
    
    Tdv = Tdv.sort_values('DiffDist', ascending=False).reset_index(drop=True)
    
    if fname is not None:
        filename = f"{fname}_diff_var.csv"
        Tdv.to_csv(filename, index=False)
        
    return Tdv