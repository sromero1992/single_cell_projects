# Author: Selim Romero, Texas A&M University
"""
Main pipeline orchestrating the full QUBO differential co-expression workflow.

Identifies a K-gene subnetwork within a biological pathway whose pairwise
co-expression has shifted most between two cellular conditions (condition A =
test/experimental; condition B = reference/control).

Condition labels are fully general: disease vs. healthy, knockout vs. wild
type, treated vs. untreated, one cell type vs. another, etc.
"""

import numpy as np
from qubo_dr.preprocess import (
    normalize_libsize,
    subset_by_condition,
    compute_gene_similarity,
    compute_differential,
)
from qubo_dr.graph import build_mnn_adjacency
from qubo_dr.pathway import get_pathway_genes
from qubo_dr.qubo import (
    build_pathway_network_matrix,
    assemble_qubo_matrix,
    solve_qubo_simulated_annealing,
    extract_subnetwork,
)
from qubo_dr.plot import (
    plot_coexpr_heatmap,
    plot_gene_network,
    plot_condition_heatmaps,
)


def run_pipeline(X, g, batch_id, genelist=None, K=30,
                 source=None, pathway_id=None, organism='human',
                 cond_a_label='A', cond_b_label='B',
                 method='mnn', n_neighbors=30, n_svd=50,
                 cell_state_A=None, cell_state_B=None,
                 use_pathway_prior=False,
                 penalty_scale=10, num_reads=1000,
                 plotit=True, edge_pct=95,
                 outfile='qubo_genes_solution.txt'):
    """
    End-to-end QUBO differential co-expression pipeline.

    Identifies the K-gene subnetwork within a pathway whose pairwise
    co-expression relationships have changed most from condition B
    (reference) to condition A (test).

    Parameters
    ----------
    X : np.ndarray, shape (G_all, N)
        Normalized scRNA-seq count matrix (all genes × all cells).
        Library-size normalization + log1p is recommended upstream.
    g : list of str
        Gene names, length G_all.
    batch_id : array-like of str, length N
        Condition label per cell.  Must contain entries matching
        cond_a_label and cond_b_label.
    genelist : list of str, optional
        Pathway gene list.  If None, fetched from source + pathway_id.
    K : int
        Target subnetwork size (default 30).
    source : str, optional
        Pathway database: 'kegg' or 'gobp'.  Required if genelist is None.
    pathway_id : str, optional
        KEGG pathway ID (e.g., 'hsa04310') or GO term name.
        Required if genelist is None.
    organism : str, optional
        Organism for pathway mapping (default 'human' → KEGG 'hsa').
    cond_a_label : str
        Label identifying condition A cells in batch_id (test / experimental).
        Default 'A'.
    cond_b_label : str
        Label identifying condition B cells in batch_id (reference / control).
        Default 'B'.
    method : str
        Graph construction method: 'knn' or 'mnn' (default 'mnn').
    n_neighbors : int
        Number of neighbours for KNN/MNN (default 30).
    n_svd : int
        Truncated SVD components for the gene embedding used in MNN
        (default 50).  SVD is performed on cells so each gene obtains a
        d-dimensional coordinate in cell-space; KNN is then between genes.
    cell_state_A : np.ndarray, shape (N_A,), optional
        Per-cell biological state scalar for condition A cells.
        Can be any continuous per-cell measure: UCell pathway activity
        score, pseudotime coordinate, cell potency score, differentiation
        rank, etc.  Pass None to omit the V_diff linear bias term.
    cell_state_B : np.ndarray, shape (N_B,), optional
        Same as cell_state_A but for condition B cells.
    use_pathway_prior : bool
        Whether to include the pathway membership prior matrix X_net.
        Set True only when extra candidate genes (outside the core pathway)
        are appended to genelist.  When X is already restricted to pathway
        genes this prior is a uniform constant and has no effect; it is
        therefore False by default.
    penalty_scale : float
        QUBO cardinality penalty multiplier: P = penalty_scale * max|Q1|.
        Default 10.
    num_reads : int
        Number of simulated annealing iterations (default 1000).
    plotit : bool
        Generate diagnostic plots (default True).
    edge_pct : float
        Percentile threshold for network edge display (default 95).
    outfile : str
        Output filename for selected genes (default 'qubo_genes_solution.txt').

    Returns
    -------
    results : dict
        'sub_g_net'         – list of selected gene names (length K)
        'sub_Q_net'         – subnetwork differential co-expression (K×K)
        'sub_Qv'            – per-gene contribution scores (K,)
        'selected_idx'      – boolean mask of selected genes (G,)
        'G_graph'           – networkx.Graph or None
        'Xdiff'             – differential co-expression matrix (G×G)
        'Q1'                – unpenalised QUBO matrix (G×G)
        'Q'                 – penalised QUBO matrix (G×G)
        'gene_names_subset' – all pathway-filtered gene names
        'best_energy'       – best QUBO objective value
        'figures'           – list of (fig, ax) if plotit=True
    """
    X = np.asarray(X)
    batch_id = np.asarray(batch_id, dtype=str)

    # ------------------------------------------------------------------
    # Step 1 – Pathway gene set
    # ------------------------------------------------------------------
    if genelist is None:
        if source is None or pathway_id is None:
            raise ValueError(
                "Provide genelist OR both source and pathway_id. "
                "source: 'kegg' or 'gobp'; pathway_id: KEGG ID or GO term."
            )
        print(f"Fetching {source.upper()} pathway: {pathway_id}")
        genelist = get_pathway_genes(source=source,
                                     pathway_id_or_name=pathway_id,
                                     organism=organism)
        if not genelist:
            raise ValueError(f"No genes returned for pathway '{pathway_id}'")
        print(f"  {len(genelist)} genes retrieved")

    g_upper = np.array([gene.upper() for gene in g])
    genelist_upper = [gene.upper() for gene in genelist]

    # ------------------------------------------------------------------
    # Step 2 – Restrict X to pathway genes  →  X ∈ R^(G × N)
    # ------------------------------------------------------------------
    mask = np.isin(g_upper, genelist_upper)
    X_path = X[mask, :]
    g_path = [g[i] for i in np.where(mask)[0]]
    G = X_path.shape[0]
    print(f"Pathway genes in data: {G}/{len(genelist)}")
    if G == 0:
        raise ValueError("No pathway genes found in gene list g.")

    # ------------------------------------------------------------------
    # Step 3 – Partition by condition  →  X_A ∈ R^(G × N_A), X_B ∈ R^(G × N_B)
    # ------------------------------------------------------------------
    print(f"Partitioning: '{cond_a_label}' (test) vs '{cond_b_label}' (reference)")
    X_A, idx_A = subset_by_condition(X_path, np.arange(X_path.shape[1]),
                                      batch_id, cond_a_label)
    X_B, idx_B = subset_by_condition(X_path, np.arange(X_path.shape[1]),
                                      batch_id, cond_b_label)
    print(f"  Condition A: {X_A.shape[1]} cells | Condition B: {X_B.shape[1]} cells")

    # ------------------------------------------------------------------
    # Step 4 – Non-negative cosine similarity matrices (Gram matrices)
    # ------------------------------------------------------------------
    print("Computing cosine similarity matrices...")
    SA, XA_norm = compute_gene_similarity(X_A)   # S_A = X_A_norm · X_A_norm^T
    SB, XB_norm = compute_gene_similarity(X_B)   # S_B = X_B_norm · X_B_norm^T

    # ------------------------------------------------------------------
    # Step 5 – Differential co-expression and V_diff
    # ------------------------------------------------------------------
    print("Computing differential co-expression (X_diff = S_B - S_A)...")
    cs_A = cell_state_A[idx_A] if cell_state_A is not None else None
    cs_B = cell_state_B[idx_B] if cell_state_B is not None else None

    Xdiff, Vdiff = compute_differential(
        SB, SA,
        Xwt_norm=XB_norm, Xko_norm=XA_norm,
        cs_wt=cs_B, cs_ko=cs_A,
    )

    if cell_state_A is None:
        print("  V_diff term: disabled (no cell-state scalar provided)")
    else:
        print("  V_diff term: enabled (cell-state scalar projected per gene)")

    # ------------------------------------------------------------------
    # Step 6 – MNN adjacency matrices  (gene-space, via cell-SVD embedding)
    # ------------------------------------------------------------------
    print(f"Building {method.upper()} adjacency matrices (SVD d={n_svd})...")
    Xmnn_A = build_mnn_adjacency(X_A, method=method, K=n_neighbors, n_svd=n_svd)
    Xmnn_B = build_mnn_adjacency(X_B, method=method, K=n_neighbors, n_svd=n_svd)

    # ------------------------------------------------------------------
    # Step 7 – Optional pathway prior  X_net
    # ------------------------------------------------------------------
    Xnet = None
    if use_pathway_prior:
        print("Building pathway membership prior (X_net)...")
        Xnet = build_pathway_network_matrix(g_path, genelist_upper)
    else:
        print("Pathway prior X_net: skipped (pure pathway-gene analysis)")

    # ------------------------------------------------------------------
    # Step 8 – Assemble Q1 and Q
    # ------------------------------------------------------------------
    print("Assembling QUBO matrix...")
    Q, Q1, P = assemble_qubo_matrix(
        Xdiff, Vdiff, Xmnn_A, Xmnn_B, K,
        Xnet_target=Xnet,
        penalty_scale=penalty_scale,
    )

    # ------------------------------------------------------------------
    # Step 9 – Solve
    # ------------------------------------------------------------------
    print(f"Solving QUBO (simulated annealing, {num_reads} reads)...")
    selected_idx, best_energy, _ = solve_qubo_simulated_annealing(
        Q, num_reads=num_reads, seed=42)

    n_sel = selected_idx.sum()
    print(f"  Selected {n_sel} genes (target K={K}) | energy={best_energy:.4f}")

    # ------------------------------------------------------------------
    # Step 10 – Extract subnetwork
    # ------------------------------------------------------------------
    sub_Q_net, sub_Qv = extract_subnetwork(Xdiff, Q1, selected_idx)
    sub_g_net = [g_path[i] for i in np.where(selected_idx)[0]]

    # ------------------------------------------------------------------
    # Step 11 – Save
    # ------------------------------------------------------------------
    _save_results(sub_g_net, outfile)

    # ------------------------------------------------------------------
    # Step 12 – Plots
    # ------------------------------------------------------------------
    figures = []
    G_graph = None
    if plotit:
        print("Generating plots...")
        figures.append((plot_condition_heatmaps(SA, SB, Xdiff, g_path), None))
        fig_s, ax_s = plot_coexpr_heatmap(
            sub_Q_net, sub_g_net,
            title=f'Subnetwork Co-expression (K={n_sel})')
        figures.append((fig_s, ax_s))
        fig_n, ax_n, G_graph = plot_gene_network(
            sub_Q_net, sub_g_net, edge_pct=edge_pct,
            title=f'Gene Network ({n_sel} genes)')
        figures.append((fig_n, ax_n))

    print("Pipeline complete.")
    return {
        'sub_g_net':          sub_g_net,
        'sub_Q_net':          sub_Q_net,
        'sub_Qv':             sub_Qv,
        'selected_idx':       selected_idx,
        'G_graph':            G_graph,
        'Xdiff':              Xdiff,
        'Q1':                 Q1,
        'Q':                  Q,
        'gene_names_subset':  g_path,
        'best_energy':        best_energy,
        'figures':            figures,
        'n_selected':         n_sel,
        'n_pathway_genes':    len(genelist),
    }


def _save_results(gene_names, outfile):
    with open(outfile, 'w') as f:
        f.write("# QUBO differential co-expression — selected subnetwork\n")
        f.write(f"# N selected genes: {len(gene_names)}\n#\n")
        for gene in gene_names:
            f.write(f"{gene}\n")
    print(f"  Saved {len(gene_names)} genes to {outfile}")
