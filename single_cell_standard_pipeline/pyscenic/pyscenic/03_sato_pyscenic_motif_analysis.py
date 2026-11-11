#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import matplotlib.pyplot as plt


# In[2]:


wdir = r"Z:\selim_working_dir\2025_sato_anestacia_circadian_rhythm\r_pre_process\SCENIC_Output"
os.chdir( wdir )
dbs = r"\cisTarget_databases"
fname = r"sato_mouse_breast_cancer_scevan_decontX_QC_annotated_10_16_2025_subann.h5ad"
fname_metadata = r"metadata_circadian_mouse_breast_decontX_QC_10_15_2025_subann_matlab.xlsx"
fname_pyscenic_output = r"pyscenic_output.loom"


# In[3]:


adata = sc.read_h5ad(fname)
df = pd.read_excel(fname_metadata)


# In[4]:


# This is necessary for a direct comparison with adata.obs.index
df = df.set_index('CellID')
assert len(adata.obs) == len(df), "The number of cells in adata and df are different."
assert all(adata.obs.index == df.index), "Cell IDs in adata.obs and df are not identical or are not in the same order."
print("Assertion successful: Cell IDs in adata.obs and df match perfectly.")
# Select all columns from df *except* the index (which is 'CellID')
new_metadata_cols = df.columns

# Use the .loc indexer to assign the columns directly
# Since the indices already match (as asserted above), this is safe and fast.
adata.obs.loc[:, new_metadata_cols] = df

# Backup the raw count matrix...
adata.layers['counts'] = adata.X



# In[5]:


from scipy.sparse import issparse
def save_loom_file( adata, fname, layer: str = None):
    """
    Saves an anndata object to a loom file.
    
    WARNING: The loom file will be saved in a non-standard (Genes x Cells) 
    orientation due to the required matrix transposition.
    """
    # Select the matrix to be saved (X or a layer)
    if layer is None:
        mat = adata.X
    else:
        # Correct access to the layer matrix
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers.")
        mat = adata.layers[layer]

    # --- Row Attributes (Genes) ---
    row_attrs = { 
        "Gene": np.array(adata.var.index) ,
    }

    if issparse(mat): # Use the imported function to check the type
        mat = mat.toarray()

    # --- Column Attributes (Cells) ---
    # NOTE: Attribute calculation must be done on the standard (Cells x Genes) matrix.
    # The result (N_cells length) is then assigned to col_attrs.
    # Since the final loom matrix will be transposed (Genes x Cells), 
    # the col_attrs will be incorrectly associated with the gene axis, 
    # but the calculation itself is correct for the cell metrics.
    
    # Calculate cell metrics on the selected matrix 'mat'
    # Use axis=1 to sum across genes for each cell
    n_gene_counts = np.array((mat > 0).sum(axis=1)).flatten()
    n_umi_counts = np.array(mat.sum(axis=1)).flatten()

    col_attrs = { 
        "CellID": np.array(adata.obs.index),
        "nGene": n_gene_counts,
        "nUMI": n_umi_counts,
    }

    # CRITICAL WARNING: Retaining the non-standard transposition as requested.
    # The matrix passed to lp.create is (Genes x Cells).
    lp.create(fname, mat.transpose(), row_attrs, col_attrs)

    return


# In[6]:


new_order = ['Early', 'Intermediate', 'Advanced']
adata.obs['tumor_stage'] = adata.obs['tumor_stage'].cat.reorder_categories(
    new_order, ordered=True)
ordered_categories = adata.obs['tumor_stage'].cat.categories
cancer_stages = ordered_categories.tolist()
print(cancer_stages)

cell_levels = [
    "Adipocytes", "B cells", "Basal KC-like", "Endothelial",
    "Epithelial cells", "Fibroblasts", "Macrophages", "Mast cells",
    "NK cells", "Neutrophils", "T cells", "Tumor", "pDCs", "Myocytes"
]
adata.obs['broad_cell_types'] = adata.obs['broad_cell_types'].cat.reorder_categories(
    cell_levels, ordered=True)
ordered_categories = adata.obs['broad_cell_types'].cat.categories
cell_types = ordered_categories.tolist()
print(cell_types)


# In[ ]:


import os
import glob
import numpy as np # Ensure numpy is imported for the filtering section

# --- Configuration (FIXED PATHS) ---
dbs_dir = "cisTarget_databases" 
f_tfs = os.path.join(dbs_dir, 'allTFs_mm.txt')
nprocs = 10

# Use os.path.join for robustness
f_db_glob = os.path.join(dbs_dir, "*.feather")
f_db_names = ' '.join( glob.glob(f_db_glob) )
f_motif_path = os.path.join(dbs_dir, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
# ------------------------------------
for icells_type in cell_types:
    icells = adata.obs['broad_cell_types'] == icells_type
    icells_type_sanitized = icells_type.replace(" ", "_")
    # CRITICAL CHECK: Add check for zero total cells for the broad type
    if sum(icells) == 0:
        print(f"Not enough cells...")
        continue
    
    print(f"\n--- Cell Type: {icells_type} (Total Cells: {sum(icells)}) ---")
    
    for istage in cancer_stages:
        ican_stage = adata.obs['tumor_stage'] == istage
        
        idx = ican_stage & icells
        num_cells = sum(idx) # Get current cell count
        
        # Check if the subset has cells before proceeding
        if num_cells == 0:
            print(f"  Stage '{istage}': 0 cells. Skipping.")
            continue
            
        sub_adata = adata[idx, :].copy()
        
        # --- GENE FILTERING ---
        min_cells_percent = 0.10
        min_cells = int(num_cells * min_cells_percent) 

        mat_to_filter = sub_adata.layers['counts'] if 'counts' in sub_adata.layers else sub_adata.X
        gene_counts = (mat_to_filter > 0).sum(axis=0).A.flatten() 
        genes_to_keep = gene_counts >= min_cells
        
        # Filter the adata object
        sub_adata = sub_adata[:, genes_to_keep].copy()
        n_genes = sub_adata.shape[1]
        
        print(f"  Stage '{istage}': {num_cells} cells and {n_genes} genes (Filtered)")

        # --- PYSCENIC STEPS ---
        fname_loom = f"{icells_type_sanitized}_{istage}.loom"

        # 1. Save Loom File
        save_loom_file( sub_adata, fname_loom, layer='counts')

        # 2. PySCENIC GRN
        fname_adj = f"{icells_type_sanitized}_{istage}_adj.csv"
        print("  Running GRN...")
        get_ipython().system('pyscenic grn {fname_loom} {f_tfs} -o {fname_adj} --num_workers {nprocs}')
        
        # 3. PySCENIC CTX (FIXED OUTPUT PARAMETER)
        fname_reg = f"{icells_type_sanitized}_{istage}_regulons.csv" # Renamed output for clarity
        print("  Running CTX...")
        get_ipython().system('pyscenic ctx {fname_adj} {f_db_names} --annotations_fname {f_motif_path} --expression_mtx_fname {fname_loom} --output {fname_reg} --mask_dropouts --num_workers {nprocs}')

