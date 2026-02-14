import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.stats import ttest_ind
import scanpy as sc
import os

from sc_splinedv import sc_splinedv, sc_norm

# ==========================================
#  USER CONFIGURATION
# ==========================================

# 1. Path to your .h5ad file
FILE_PATH = "data/GSM6061702_GSM6061683.h5ad"  

# 2. Output folder
OUTPUT_FOLDER = "results/"

# 3. Column names in your H5AD
COLUMN_CELLTYPE = "CellType" 
COLUMN_BATCH = "BatchID" 

# 4. Groups
# Batch 1 = Control (Reference), Batch 2 = Cancer (Target)
BATCH_1_NAME = "Unaffected"   
BATCH_2_NAME = "CRC"          

# 5. Target Cell Type
TARGET_CELL_TYPE = 'Epithelial cells'
# ==========================================

def sc_filter_cells_only(X, genes, max_mito_pct=0.05):
    """Filters CELLS based on mitochondrial content. Keeps all genes."""
    if sparse.issparse(X):
        X = X.toarray()
    
    genes = np.array(genes)
    is_mito = np.array([str(g).upper().startswith('MT-') for g in genes])
    
    if np.any(is_mito):
        mito_counts = np.sum(X[is_mito, :], axis=0)
        total_counts = np.sum(X, axis=0)
        mito_pct = mito_counts / (total_counts + 1e-10)
        cell_mask = mito_pct <= max_mito_pct
        X = X[:, cell_mask]
    
    return X

def bh_adjust(pvals):
    """Calculates Benjamini-Hochberg FDR correction manually."""
    pvals = np.array(pvals)
    n = len(pvals)
    
    # Sort p-values
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    
    # Calculate FDR: p * N / rank
    ranks = np.arange(1, n + 1)
    adjusted = sorted_pvals * n / ranks
    
    # 1. Cap at 1.0
    adjusted = np.minimum(adjusted, 1.0)
    
    # 2. Enforce monotonicity (p_adj[i] <= p_adj[i+1])
    # Iterate backwards and keep the minimum
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i+1])
        
    # Map back to original order
    results = np.zeros(n)
    results[sorted_indices] = adjusted
    return results

def sc_deg(X1, X2, genes, name1='Set1', name2='Set2'):
    """ 
    Perform Differential Expression (Welch's T-test) + FDR Correction.
    Matches Set A column format exactly.
    """
    genes = np.array(genes)
    n_genes = len(genes)
    
    # 1. Calculate Means
    mean1 = np.mean(X1, axis=1)
    mean2 = np.mean(X2, axis=1)

    # 2. Calculate Percentages
    pct1 = np.sum(X1 > 0, axis=1) / X1.shape[1]
    pct2 = np.sum(X2 > 0, axis=1) / X2.shape[1]
    
    # 3. Log Fold Change (log2(Target / Reference))
    log_fc = np.log2((mean2 + 0.1) / (mean1 + 0.1))
    
    # 4. T-test
    pvals = np.zeros(n_genes)
    for i in range(n_genes):
        try:
            if np.var(X1[i, :]) == 0 and np.var(X2[i, :]) == 0:
                pvals[i] = 1.0
            else:
                _, pval = ttest_ind(X1[i, :], X2[i, :], equal_var=False)
                # Handle NaNs or underflow
                pvals[i] = pval if not np.isnan(pval) else 1.0
        except:
            pvals[i] = 1.0
            
    # 5. FDR Correction (THE FIX)
    p_adj = bh_adjust(pvals)

    # Create Table with Set A style headers
    Tde = pd.DataFrame({
        'gene': genes,
        'p_val': pvals,          
        'avg_log2FC': log_fc,
        'abs_log2FC': np.abs(log_fc),
        f'avg_1_{name2}': mean2,       
        f'avg_2_{name1}': mean1,
        f'pct_1_{name2}': pct2,
        f'pct_2_{name1}': pct1,
        'p_val_adj': p_adj  # Matches Set A
    })
    
    # Sort by p-value
    Tde = Tde.sort_values('p_val').reset_index(drop=True)
    
    return Tde

def main():
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    print(f"Loading data from {FILE_PATH}...")
    try:
        adata = sc.read_h5ad(FILE_PATH)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    if 'counts' in adata.layers:
        X_all = adata.layers['counts']
    elif adata.raw is not None:
        X_all = adata.raw.X
    else:
        X_all = adata.X

    if adata.raw is not None:
         genes_all = adata.raw.var_names.values
    else:
         genes_all = adata.var_names.values

    cell_types = adata.obs[COLUMN_CELLTYPE].values.astype(str)
    batch_ids = adata.obs[COLUMN_BATCH].values.astype(str)

    print(f"Analyzing: {TARGET_CELL_TYPE} ({BATCH_1_NAME} vs {BATCH_2_NAME})")
    
    # --- FILTERING (UNION METHOD) ---
    print(f"Processing {BATCH_1_NAME}...")
    idx1 = (batch_ids == BATCH_1_NAME) & (cell_types == TARGET_CELL_TYPE)
    if np.sum(idx1) == 0: print("Error: No cells for Batch 1"); return
    X1_clean = sc_filter_cells_only(X_all[idx1, :].T, genes_all)
    
    print(f"Processing {BATCH_2_NAME}...")
    idx2 = (batch_ids == BATCH_2_NAME) & (cell_types == TARGET_CELL_TYPE)
    if np.sum(idx2) == 0: print("Error: No cells for Batch 2"); return
    X2_clean = sc_filter_cells_only(X_all[idx2, :].T, genes_all)
    
    print("Selecting genes (Union)...")
    gene_mask = (np.sum(X1_clean > 0, axis=1) >= 1) | (np.sum(X2_clean > 0, axis=1) >= 1)
    
    gl = genes_all[gene_mask]
    X1 = X1_clean[gene_mask, :]
    X2 = X2_clean[gene_mask, :]
    print(f"Genes: {len(gl)}")
    
    # --- ANALYSIS ---
    label = f"{TARGET_CELL_TYPE.replace(' ', '_')}-{BATCH_1_NAME}_vs_{BATCH_2_NAME}"
    
    print("Running Analysis...")
    X1_norm = sc_norm(X1, 'libsize')
    X2_norm = sc_norm(X2, 'libsize')
    
    # Run DE
    Tde = sc_deg(X1_norm, X2_norm, gl, name1=BATCH_1_NAME, name2=BATCH_2_NAME)
    
    # Run DV (Optional, reusing norm data)
    fname_base = os.path.join(OUTPUT_FOLDER, label)
    Tdv = sc_splinedv(X1_norm, X2_norm, gl, gl, 
                      name1=BATCH_1_NAME, name2=BATCH_2_NAME, 
                      fname=fname_base, plotit=False)
    
    # Save Results
    Tup = Tde[(Tde['p_val'] < 0.05) & (Tde['avg_log2FC'] > 0.5)]
    Tdn = Tde[(Tde['p_val'] < 0.05) & (Tde['avg_log2FC'] < -0.5)]
    
    fname_de = os.path.join(OUTPUT_FOLDER, f"DE_{label}.xlsx")
    print(f"Saving to {fname_de}...")
    
    with pd.ExcelWriter(fname_de, engine='openpyxl') as writer:
        Tde.to_excel(writer, sheet_name='All genes', index=False)
        Tup.to_excel(writer, sheet_name='Up-regulated', index=False)
        Tdn.to_excel(writer, sheet_name='Down-regulated', index=False)
        
    print(f"Completed. Top DE Gene: {Tde.iloc[0]['gene']}, P-adj: {Tde.iloc[0]['p_val_adj']}")

if __name__ == "__main__":
    main()