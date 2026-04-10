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
    compute_spectral_penalty,
    solve_qubo_simulated_annealing,
    extract_subnetwork,
    solve_classical_topk,
)
from qubo_dr.plot import (
    plot_coexpr_heatmap,
    plot_gene_network,
    plot_condition_heatmaps,
)
from qubo_dr.permutation import test_differential_hub


def run_pipeline(X, g, batch_id,
                 alternative_pathway_genelist=None,
                 K=30,
                 pathway_name=None, library=None, organism='human',
                 # backward-compat aliases
                 source=None, pathway_id=None,
                 genelist=None,          # deprecated: use alternative_pathway_genelist
                 cond_a_label=None, cond_b_label=None,
                 test_label='disease', reference_label='control',
                 method='mnn', n_neighbors=30, n_svd=50,
                 extra_genes=None,
                 forced_genes=None,
                 cell_state=None,
                 use_pathway_prior=None,
                 penalty_scale=2, num_reads=1000,
                 solver='sa',
                 plotit=True,
                 top_pct=10,
                 run_perm_test=True, n_perm=200,
                 outfile='qubo_genes_solution.txt'):
    """
    End-to-end QUBO differential co-expression pipeline.

    Finds the K-gene hub inside a biological pathway with maximum pairwise
    co-expression shift from reference (B) → test (A) via QUBO optimisation.

    Parameters
    ----------
    X : ndarray (G_all, N)  — normalised scRNA-seq matrix (genes × cells).
    g : list[str]           — gene names, length G_all.
    batch_id : array[str]   — condition label per cell (length N).
    alternative_pathway_genelist : list[str] or None
        Explicit gene list. If None, fetched from Enrichr via pathway_name.
    K : int                 — hub size (default 30).
    pathway_name : str      — Enrichr keyword, e.g. 'Wnt signaling pathway'.
    library : str or None   — restrict Enrichr lookup to one library.
    organism : str          — 'human' or 'mouse' (default 'human').
    cond_a_label : str      — condition A label in batch_id (test, default 'A').
    cond_b_label : str      — condition B label in batch_id (ref,  default 'B').
    test_label : str        — display name for condition A (default 'disease').
    reference_label : str   — display name for condition B (default 'control').
    method : str            — neighbor graph: 'mnn' (default) or 'knn'.
    n_neighbors : int       — K for MNN/KNN graph (default 30).
    n_svd : int             — SVD dims for gene embedding before MNN (default 50).
    cell_state : ndarray (N,) or None
        Per-cell biological state scalar (UCell score, pseudotime, potency,
        etc.) for all cells — split by condition internally.
        Enables the V_diff linear bias term. Pass None to omit (default).
    extra_genes : list[str] or None
        Genes added to the pathway pool to compete for K slots.
        Enables pathway-membership prior (Xnet) automatically.
    forced_genes : list[str] or None
        Genes guaranteed in the hub regardless of optimisation.
        Implementation: variable reduction — forces z_i=1, solves only the
        remaining K−|forced| free slots.  Missing genes are skipped with a
        warning.  Result is always exactly K genes.
    use_pathway_prior : bool or None
        Include Xnet (pathway-membership prior). None = auto (True when
        extra_genes provided, False for pure pathway analysis).
    penalty_scale : float   — spectral penalty multiplier P = scale·‖Q1‖₂ (default 2).
    num_reads : int         — SA restarts (default 1000).
    solver : str            — 'sa' (neal), 'sd' (steepest-descent), 'fallback' (NumPy SA).
    plotit : bool           — generate diagnostic plots (default True).
    top_pct : float         — % of edges shown in Q0 network plot (default 10).
    run_perm_test : bool    — hub permutation test (default True).
    n_perm : int            — permutations for hub test (default 200).
    outfile : str           — output filename for selected genes.

    Returns
    -------
    dict:  sub_g_net, sub_Q0, sub_Q_net, sub_Q1, sub_Qv, selected_idx,
           dS, Q0, Q1, Q, best_energy, hub_perm_result, gene_names_subset,
           forced_genes_used, figures, n_selected, n_pathway_genes.
    """
    X = np.asarray(X)
    batch_id = np.asarray(batch_id, dtype=str)

    # ── Backward compatibility ─────────────────────────────────────────────
    import warnings
    if genelist is not None and alternative_pathway_genelist is None:
        warnings.warn(
            "The 'genelist' parameter is deprecated. "
            "Use 'alternative_pathway_genelist' instead.",
            DeprecationWarning, stacklevel=2,
        )
        alternative_pathway_genelist = genelist
    if cond_a_label is not None:
        test_label = cond_a_label
    if cond_b_label is not None:
        reference_label = cond_b_label
    if pathway_id is not None and pathway_name is None:
        pathway_name = pathway_id
    if source is not None:
        warnings.warn(
            "The 'source' argument is deprecated and ignored. "
            "DR-QUBO now uses Enrichr (via gseapy) for all pathway lookups. "
            "Pass 'library' to restrict to a specific Enrichr library.",
            DeprecationWarning, stacklevel=2,
        )

    # ------------------------------------------------------------------
    # Step 1 – Pathway gene set
    # ------------------------------------------------------------------
    if alternative_pathway_genelist is None:
        if pathway_name is None:
            raise ValueError(
                "Provide alternative_pathway_genelist OR pathway_name.\n"
                "  pathway_name               : keyword or exact term, e.g. 'Wnt signaling pathway'\n"
                "  library                    : optional Enrichr library, e.g. 'KEGG_2026'\n"
                "                               (if omitted, all default libraries are searched)\n"
                "  alternative_pathway_genelist: explicit gene list — skips Enrichr lookup\n"
                "  Tip: call browse_pathways('Wnt') first to browse available terms."
            )
        lib_tag = f" [{library}]" if library else " [auto-search: KEGG, GO, Reactome, ...]"
        print(f"Fetching pathway '{pathway_name}'{lib_tag}")
        alternative_pathway_genelist = get_pathway_genes(
            pathway_name=pathway_name, library=library, organism=organism)
        if not alternative_pathway_genelist:
            raise ValueError(f"No genes returned for pathway '{pathway_name}'")
        print(f"  {len(alternative_pathway_genelist)} genes retrieved")

    # Core pathway gene set (uppercase, used for prior matrix later)
    genelist_upper = [gene.upper() for gene in alternative_pathway_genelist]

    # ── Append extra genes (GWAS hits, drug targets, TFs, …) ──────────────
    if extra_genes:
        pathway_set = set(genelist_upper)
        new_extras  = [eg for eg in extra_genes
                       if eg.upper() not in pathway_set]
        if new_extras:
            print(f"  Appending {len(new_extras)} extra genes: {new_extras}")
        genelist_full = list(alternative_pathway_genelist) + new_extras
    else:
        genelist_full = list(alternative_pathway_genelist)

    genelist_full_upper = [gene.upper() for gene in genelist_full]

    # ── Auto-enable pathway prior when extras are present ─────────────────
    # Pure pathway:  Xnet is -1 uniformly off-diagonal = constant energy shift
    #                → no effect on gene selection, skip it.
    # Mixed pathway + extras: Xnet assigns -1 bonus to pathway-pathway pairs
    #                only, steering the solver toward core pathway genes.
    #                → must be enabled.
    if use_pathway_prior is None:
        use_pathway_prior = bool(extra_genes)
    new_extras = new_extras if extra_genes else []

    g_upper = np.array([gene.upper() for gene in g])

    # ------------------------------------------------------------------
    # Step 2 – Restrict X to pathway + extra genes  →  X ∈ R^(G × N)
    #           Works for any X size (500, 20k, etc.) via boolean mask.
    # ------------------------------------------------------------------
    mask   = np.isin(g_upper, genelist_full_upper)
    X_path = X[mask, :]
    g_path = [g[i] for i in np.where(mask)[0]]
    G      = X_path.shape[0]

    n_pathway_found = sum(1 for gn in g_path if gn.upper() in set(genelist_upper))
    n_extra_found   = G - n_pathway_found
    print(f"Pathway genes in data: {n_pathway_found}/{len(alternative_pathway_genelist)}"
          + (f"  |  extra genes: {n_extra_found}/{len(new_extras) if extra_genes else 0}"
             if extra_genes else ""))
    if G == 0:
        raise ValueError("No pathway genes found in gene list g.")

    # ------------------------------------------------------------------
    # Step 3 – Partition by condition  →  X_A ∈ R^(G × N_A), X_B ∈ R^(G × N_B)
    # ------------------------------------------------------------------
    print(f"Partitioning: '{test_label}' (test/A) vs '{reference_label}' (reference/B)")
    X_A, idx_A = subset_by_condition(X_path, np.arange(X_path.shape[1]),
                                      batch_id, test_label)
    X_B, idx_B = subset_by_condition(X_path, np.arange(X_path.shape[1]),
                                      batch_id, reference_label)
    print(f"  Test (A): {X_A.shape[1]} cells | Reference (B): {X_B.shape[1]} cells")

    # ------------------------------------------------------------------
    # Step 4 – Non-negative cosine similarity matrices (Gram matrices)
    # ------------------------------------------------------------------
    print("Computing cosine similarity matrices...")
    SA, XA_norm = compute_gene_similarity(X_A)   # S_A = X_A_norm · X_A_norm^T
    SB, XB_norm = compute_gene_similarity(X_B)   # S_B = X_B_norm · X_B_norm^T

    # ------------------------------------------------------------------
    # Step 5 – Differential co-expression and V_diff
    # ------------------------------------------------------------------
    print("Computing differential co-expression (dS = S_ref − S_test)...")
    cs_A = cell_state[idx_A] if cell_state is not None else None
    cs_B = cell_state[idx_B] if cell_state is not None else None

    dS, Vdiff = compute_differential(
        SB, SA,
        Xref_norm=XB_norm, Xtest_norm=XA_norm,
        cs_ref=cs_B, cs_test=cs_A,
    )

    if cell_state is None:
        print("  V_diff term: disabled (no cell-state scalar provided)")
    else:
        print("  V_diff term: enabled (cell-state scalar projected per gene)")

    # ------------------------------------------------------------------
    # Step 6 – MNN adjacency matrices  (gene-space, via cell-SVD embedding)
    # ------------------------------------------------------------------
    print(f"Building {method.upper()} adjacency matrices (SVD d={n_svd})...")
    A_mnn_test = build_mnn_adjacency(X_A, method=method, K=n_neighbors, n_svd=n_svd)
    A_mnn_ref = build_mnn_adjacency(X_B, method=method, K=n_neighbors, n_svd=n_svd)

    # ------------------------------------------------------------------
    # Step 7 – Optional pathway prior  X_net
    # ------------------------------------------------------------------
    Xnet = None
    if use_pathway_prior:
        n_extras = len(genelist_full) - len(alternative_pathway_genelist)
        print(f"Building pathway membership prior (X_net)  "
              f"[{len(alternative_pathway_genelist)} pathway + {n_extras} extra genes]...")
        # genelist_upper = core pathway only → assigns -1 to pathway-pathway
        # pairs, 0 to extra-* pairs, giving the solver a bias toward selecting
        # core pathway genes when extras compete for the K slots.
        Xnet = build_pathway_network_matrix(g_path, genelist_upper)
    else:
        print("Pathway prior X_net: skipped (pure pathway-gene analysis — "
              "prior would be a uniform constant with no effect)")

    # ------------------------------------------------------------------
    # Step 8 – Assemble Q1 and Q
    # ------------------------------------------------------------------
    print("Assembling QUBO matrix...")
    Q, Q1, Q0, P = assemble_qubo_matrix(
        dS, Vdiff, A_mnn_test, A_mnn_ref, K,
        Xnet_target=Xnet,
        penalty_scale=penalty_scale,
    )

    # ------------------------------------------------------------------
    # Step 9 – Forced-gene variable reduction (if forced_genes provided)
    # ------------------------------------------------------------------
    g_path_upper = np.array([gn.upper() for gn in g_path])
    forced_idx   = np.array([], dtype=int)   # will hold indices into g_path

    if forced_genes:
        found, missing = [], []
        for fg in forced_genes:
            hits = np.where(g_path_upper == fg.upper())[0]
            if len(hits) > 0:
                found.append(int(hits[0]))
            else:
                missing.append(fg)
        if missing:
            print(f"  [forced_genes] WARNING: {missing} not found in pathway gene set "
                  f"— skipped.")
        if found:
            forced_idx = np.unique(np.array(found, dtype=int))
            print(f"  [forced_genes] Pinning {len(forced_idx)} gene(s) to hub: "
                  f"{[g_path[i] for i in forced_idx]}")

    # ------------------------------------------------------------------
    # Step 10 – Solve
    # ------------------------------------------------------------------
    solver_label = {'sa': 'SimulatedAnnealing (neal)',
                    'sd': 'SteepestDescent (classical)',
                    'fallback': 'NumPy SA (fallback)'}.get(solver, solver)

    if len(forced_idx) > 0:
        K_free    = K - len(forced_idx)
        free_idx  = np.array([i for i in range(G) if i not in set(forced_idx.tolist())])

        if K_free <= 0:
            # All K slots consumed by forced genes
            if K_free < 0:
                print(f"  WARNING: {len(forced_idx)} forced genes > K={K}. "
                      f"Using first {K} forced genes only.")
                forced_idx = forced_idx[:K]
            print(f"  All K={K} slots filled by forced genes — skipping optimisation.")
            selected_idx = np.zeros(G, dtype=bool)
            selected_idx[forced_idx] = True
            best_energy  = float(selected_idx.astype(float) @ Q @ selected_idx.astype(float))
        else:
            # ── Variable reduction ───────────────────────────────────────────
            # Fold cross-terms (Q1_forced→free) into the free-gene diagonal:
            #   E = const + z_U^T [Q1_UU + diag(2·Q1_FU^T·1_F)] z_U  +  penalty(K')
            Q1_UU       = Q1[np.ix_(free_idx, free_idx)]
            linear_bias = 2.0 * Q1[np.ix_(forced_idx, free_idx)].sum(axis=0)

            Q1_red = Q1_UU.copy()
            np.fill_diagonal(Q1_red, Q1_red.diagonal() + linear_bias)

            # Rebuild penalty for K_free on the reduced (n_free × n_free) matrix
            P_red, sigma_red = compute_spectral_penalty(Q1_red, penalty_scale)
            n_free = len(free_idx)
            Q_red  = Q1_red.copy()
            Q_red += P_red
            np.fill_diagonal(Q_red, Q1_red.diagonal() + P_red * (1 - 2 * K_free))
            Q_red  = (Q_red + Q_red.T) / 2.0

            print(f"Solving reduced QUBO ({solver_label}, {num_reads} reads) | "
                  f"{n_free} free genes → K'={K_free} slots ...")
            sel_free, best_energy, _ = solve_qubo_simulated_annealing(
                Q_red, num_reads=num_reads, seed=42, solver=solver)

            # Reconstruct full selected_idx
            selected_idx = np.zeros(G, dtype=bool)
            selected_idx[forced_idx]           = True
            selected_idx[free_idx[sel_free]]   = True
    else:
        print(f"Solving QUBO ({solver_label}, {num_reads} reads)...")
        selected_idx, best_energy, _ = solve_qubo_simulated_annealing(
            Q, num_reads=num_reads, seed=42, solver=solver)

    n_sel = selected_idx.sum()
    print(f"  Selected {n_sel} genes (target K={K}) | energy={best_energy:.4f}")

    # ------------------------------------------------------------------
    # Step 11 – Extract subnetwork  (K×K matrices for selected genes)
    # ------------------------------------------------------------------
    sub_Q0, sub_Q1_net, sub_Qv = extract_subnetwork(dS, Q1, selected_idx)
    # sub_Q0      = K×K dS for selected genes (pure differential signal)
    # sub_Q1_net  = K×K Q1  for selected genes  (optimization landscape)
    # sub_Q_net kept as alias for backward compat (= -sub_Q0 for the heatmap)
    sub_Q_net = -sub_Q0   # negated: positive = test gain (for heatmap display)
    sub_g_net = [g_path[i] for i in np.where(selected_idx)[0]]

    # ------------------------------------------------------------------
    # Step 12 – Save
    # ------------------------------------------------------------------
    _save_results(sub_g_net, outfile)

    # ------------------------------------------------------------------
    # Step 13 – Hub permutation test (is the selected subnetwork
    #           significantly differential?)
    # ------------------------------------------------------------------
    hub_perm_result = None
    if run_perm_test:
        print(f"Hub significance test ({n_perm} permutations)...")
        hub_perm_result = test_differential_hub(
            X_A, X_B, dS, selected_idx,
            n_perm=n_perm, seed=42, verbose=True,
        )

    # ------------------------------------------------------------------
    # Step 14 – Plots
    # ------------------------------------------------------------------
    figures  = []
    G_graph0 = None   # networkx graph for Q0 network

    pval_ann   = hub_perm_result['pval']   if hub_perm_result else None
    zscore_ann = hub_perm_result['zscore'] if hub_perm_result else None

    if plotit:
        print("Generating plots...")

        # 3-panel co-expression heatmaps (SA, SB, dS)
        figures.append((plot_condition_heatmaps(
            SA, SB, dS, g_path,
            test_label=test_label, reference_label=reference_label), None))

        # K×K subnetwork heatmap
        fig_s, ax_s = plot_coexpr_heatmap(
            sub_Q_net, sub_g_net,
            title=(f'Subnetwork Differential Co-expression '
                   f'({test_label} vs {reference_label}, K={n_sel})'))
        figures.append((fig_s, ax_s))

        # ── Hub network: only the K selected genes, ALL their edges ───────
        # sub_Q0 = K×K dS  (pure biological differential signal)
        # With K≈20 genes there are at most K*(K-1)/2 ≈ 190 edges — no
        # thinning needed.  All nodes are hub members → all drawn in orange.
        all_hub = np.ones(n_sel, dtype=bool)

        fig_n0, ax_n0, G_graph0 = plot_gene_network(
            sub_Q0, sub_g_net,
            selected_mask=all_hub,
            top_pct=100,                    # show every edge — K is small
            test_label=test_label,
            reference_label=reference_label,
            perm_pval=pval_ann, perm_zscore=zscore_ann,
            title=(f'Q0 Hub (K={n_sel})  '
                   f'[{test_label} − {reference_label}]  red=test gain, blue=ref gain'))
        figures.append((fig_n0, ax_n0))

    print("Pipeline complete.")
    return {
        'sub_g_net':          sub_g_net,
        'sub_Q0':             sub_Q0,        # K×K dS for hub genes (Q0 subnetwork)
        'sub_Q1':             sub_Q1_net,    # K×K Q1 for hub genes (kept for downstream use)
        'sub_Q_net':          sub_Q_net,     # = −sub_Q0  (for heatmap: pos = test/disease gain)
        'sub_Qv':             sub_Qv,
        'selected_idx':       selected_idx,
        'G_graph0':           G_graph0,      # networkx Graph — Q0 network
        'dS':                 dS,
        'Q0':                 Q0,            # pure dS = S_ref − S_test  (full pathway G×G)
        'Q1':                 Q1,            # Q0 + Vdiff − MNN [− prior]  (full pathway G×G)
        'Q':                  Q,             # Q1 + cardinality penalty  (full pathway G×G)
        'gene_names_subset':  g_path,
        'best_energy':        best_energy,
        'hub_perm_result':    hub_perm_result,
        'figures':            figures,
        'n_selected':         n_sel,
        'n_pathway_genes':    len(alternative_pathway_genelist),
        'forced_genes_used':  [g_path[i] for i in forced_idx] if len(forced_idx) > 0 else [],
    }


def _save_results(gene_names, outfile):
    with open(outfile, 'w') as f:
        f.write("# QUBO differential co-expression — selected subnetwork\n")
        f.write(f"# N selected genes: {len(gene_names)}\n#\n")
        for gene in gene_names:
            f.write(f"{gene}\n")
    print(f"  Saved {len(gene_names)} genes to {outfile}")


# ===========================================================================
#  Classical baseline — same preprocessing, algebraic gene selection
# ===========================================================================

def run_pipeline_classical(X, g, batch_id,
                            alternative_pathway_genelist=None,
                            K=30,
                            pathway_name=None, library=None, organism='human',
                            genelist=None,   # deprecated: use alternative_pathway_genelist
                            test_label='disease', reference_label='control',
                            method='rowsum',
                            method_mnn='mnn', n_neighbors=30, n_svd=50,
                            plotit=True,
                            top_pct=10,
                            run_perm_test=True, n_perm=200,
                            outfile='classical_genes_solution.txt'):
    """
    Classical baseline for differential co-expression gene selection.

    Runs the same preprocessing as run_pipeline (library-size normalisation,
    cosine similarity, dS = S_B − S_A, optional MNN adjacency) but selects
    genes by a simple algebraic ranking instead of QUBO optimisation.

    Use this to benchmark the QUBO solution: the genes run_pipeline picks
    but run_pipeline_classical misses are the ones where the joint
    combinatorial structure (enforced by the MNN graph and the cardinality
    constraint) adds information beyond a simple per-gene score.

    Parameters
    ----------
    X, g, batch_id, alternative_pathway_genelist, K, pathway_name, library, organism :
        Same as run_pipeline.
    test_label, reference_label : str
        Condition labels.
    method : str, {'rowsum', 'absrowsum', 'spectral'}
        Classical ranking criterion (see solve_classical_topk).
        'rowsum'    — net disease co-expression gain per gene (default).
        'absrowsum' — total change magnitude, direction-agnostic.
        'spectral'  — leading eigenvector of the differential matrix (dS).
    method_mnn : str
        Graph method for adjacency ('mnn' or 'knn').  The adjacency weights
        the scores the same way Q1 is weighted in the QUBO, so the classical
        and QUBO objectives are directly comparable.
    plotit : bool
        Generate network plot (default True).
    top_pct : float
        Top percentage of edges (by |dS|) to display (default 10).
    run_perm_test : bool
        Run hub permutation test (default True).
    n_perm : int
        Number of permutations for hub test (default 200).
    outfile : str
        Output filename (default 'classical_genes_solution.txt').

    Returns
    -------
    dict with keys matching run_pipeline output:
        sub_g_net, sub_Q_net, sub_Qv, selected_idx, scores,
        dS, gene_names_subset, figures, n_selected, n_pathway_genes.
    """
    import numpy as np
    import warnings

    X = np.asarray(X)
    batch_id = np.asarray(batch_id, dtype=str)

    # ── Backward compatibility ─────────────────────────────────────────────
    if genelist is not None and alternative_pathway_genelist is None:
        warnings.warn(
            "The 'genelist' parameter is deprecated. "
            "Use 'alternative_pathway_genelist' instead.",
            DeprecationWarning, stacklevel=2,
        )
        alternative_pathway_genelist = genelist

    # ── Step 1: Pathway gene set ───────────────────────────────────────────
    if alternative_pathway_genelist is None:
        if pathway_name is None:
            raise ValueError(
                "Provide alternative_pathway_genelist or pathway_name.")
        lib_tag = f" [{library}]" if library else " [auto]"
        print(f"Fetching pathway '{pathway_name}'{lib_tag}")
        alternative_pathway_genelist = get_pathway_genes(
            pathway_name, library=library, organism=organism)

    print(f"Classical baseline  |  method='{method}'  |  K={K}")

    # ── Step 2: Subset to pathway genes ───────────────────────────────────
    g_upper    = [gi.upper() for gi in g]
    gl_upper   = [gi.upper() for gi in alternative_pathway_genelist]
    path_idx   = [i for i, gi in enumerate(g_upper) if gi in set(gl_upper)]
    g_path     = [g[i] for i in path_idx]
    X_path     = X[path_idx, :]
    n_in = len(g_path)
    print(f"Pathway genes in data: {n_in}/{len(alternative_pathway_genelist)}")

    # ── Step 3: Partition conditions ──────────────────────────────────────
    idx_A = np.where(batch_id == test_label)[0]
    idx_B = np.where(batch_id == reference_label)[0]
    print(f"Partitioning: '{test_label}' (A) vs '{reference_label}' (B)")
    print(f"  A: {len(idx_A)} cells | B: {len(idx_B)} cells")

    XA = X_path[:, idx_A].astype(float)
    XB = X_path[:, idx_B].astype(float)

    # ── Step 4: Cosine similarity ──────────────────────────────────────────
    def _cosine(M):
        norms = np.linalg.norm(M, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Mn = M / norms
        S = Mn @ Mn.T
        return np.clip(S, 0, None)

    SA = _cosine(XA)
    SB = _cosine(XB)
    dS = SB - SA                   # positive = reference gain / test loss

    # ── Step 5: MNN adjacency (same as QUBO pipeline for fair comparison) ─
    from qubo_dr.graph import build_mnn_adjacency
    A_adj = build_mnn_adjacency(X_path, method=method_mnn, K=n_neighbors, n_svd=n_svd)

    # ── Step 6: Classical selection ────────────────────────────────────────
    selected_idx, scores = solve_classical_topk(
        dS, K, adjacency=A_adj, method=method)

    n_sel     = selected_idx.sum()
    sel_int   = np.where(selected_idx)[0]            # integer indices (safe)
    sub_g_net = [g_path[i] for i in sel_int]
    sub_Q0    = dS[np.ix_(sel_int, sel_int)]      # K×K dS
    sub_Q_net = -sub_Q0                              # negated for heatmap display
    sub_Qv    = scores[selected_idx]

    print(f"  Selected {n_sel} genes (target K={K})")
    print(f"  Genes: {sub_g_net}")

    _save_results(sub_g_net, outfile)

    # ── Step 7: Hub permutation test ────────────────────────────────────────
    hub_perm_result = None
    if run_perm_test:
        print(f"Hub significance test ({n_perm} permutations)...")
        hub_perm_result = test_differential_hub(
            XA, XB, dS, selected_idx,
            n_perm=n_perm, seed=42, verbose=True,
        )

    # ── Step 8: Plots ──────────────────────────────────────────────────────
    figures  = []
    G_graph0 = None
    if plotit:
        print("Generating plots...")
        figures.append((plot_condition_heatmaps(
            SA, SB, dS, g_path,
            test_label=test_label, reference_label=reference_label), None))
        fig_s, ax_s = plot_coexpr_heatmap(
            sub_Q_net, sub_g_net,
            title=f'Classical Differential Co-expression ({method}, K={n_sel})')
        figures.append((fig_s, ax_s))
        pval_ann   = hub_perm_result['pval']   if hub_perm_result else None
        zscore_ann = hub_perm_result['zscore'] if hub_perm_result else None
        all_hub = np.ones(n_sel, dtype=bool)
        fig_n0, ax_n0, G_graph0 = plot_gene_network(
            sub_Q0, sub_g_net,
            selected_mask=all_hub,
            top_pct=100,                      # all edges — K is small
            test_label=test_label, reference_label=reference_label,
            perm_pval=pval_ann, perm_zscore=zscore_ann,
            title=(f'Classical Hub (K={n_sel}) — {method}  '
                   f'[{test_label} − {reference_label}]  red=test gain, blue=ref gain'))
        figures.append((fig_n0, ax_n0))

    print("Classical pipeline complete.")
    return {
        'sub_g_net':         sub_g_net,
        'sub_Q_net':         sub_Q_net,
        'sub_Qv':            sub_Qv,
        'selected_idx':      selected_idx,
        'scores':            scores,
        'dS':                dS,
        'gene_names_subset': g_path,
        'G_graph0':          G_graph0,
        'hub_perm_result':   hub_perm_result,
        'figures':           figures,
        'n_selected':        n_sel,
        'n_pathway_genes':   len(alternative_pathway_genelist),
    }
