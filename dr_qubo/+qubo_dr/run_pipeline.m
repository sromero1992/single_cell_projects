function results = run_pipeline(X, g, batch_id, genelist, K, opts)
% RUN_PIPELINE  QUBO-based differential co-expression subnetwork pipeline.
%
%   results = RUN_PIPELINE(X, g, batch_id, genelist, K)
%   results = RUN_PIPELINE(X, g, batch_id, genelist, K, opts)
%
% DESCRIPTION:
%   Identifies the K-gene subnetwork within a pathway whose pairwise
%   co-expression relationships have shifted most between two cellular
%   conditions (condition A = test/experimental; condition B =
%   reference/control).
%
%   The method is general: condition labels, pathways, and cell-state
%   scalars are all user-defined.  Typical contrasts include disease vs.
%   healthy, knockout vs. wild type, treated vs. untreated, or any two
%   annotated cell populations.
%
% INPUTS (required):
%   X        - Normalised count matrix, G_all × N (genes × cells).
%              Library-size normalisation + log1p recommended upstream.
%   g        - Gene name string array, length G_all.
%   batch_id - Condition label string array, length N.
%              Must contain entries matching cond_a_label and cond_b_label.
%   genelist - Pathway gene string array (core pathway set).
%   K        - Target subnetwork size (integer).
%
% INPUTS (optional struct fields in opts):
%   cond_a_label    - Label for condition A cells (default "A")
%   cond_b_label    - Label for condition B cells (default "B")
%   method          - Adjacency method: 'mnn' or 'knn' (default 'mnn')
%   n_neighbors     - Neighbours for MNN/KNN (default 15)
%   n_svd           - SVD components for gene embedding (default 50).
%                     SVD is run on cells so each gene gets coordinates
%                     in cell-space; KNN is then between genes.
%   Vdiff           - Pre-computed cell-state difference vector, G×1.
%                     If provided, skips internal computation.
%                     Can encode any per-cell scalar: UCell pathway
%                     activity, pseudotime, cell potency, etc.
%   cell_state_A    - Per-cell scalar vector for condition A, N_A×1.
%                     Used to compute V_diff internally if Vdiff not given.
%   cell_state_B    - Per-cell scalar vector for condition B, N_B×1.
%   use_pathway_prior - Include X_net pathway membership prior (default false).
%                     Set true only when extra candidate genes (outside
%                     the core pathway) are appended to genelist.  When
%                     X is restricted to pathway genes only, X_net is a
%                     uniform constant and has no effect on the solution.
%   penalty_scale   - Cardinality penalty multiplier (default 10).
%   plotit          - Generate diagnostic plots (default true).
%   outfile         - Output filename (default 'qubo_genes_sol.txt').
%   edge_pct        - Edge percentile threshold for network plot (default 95).
%
% OUTPUT:
%   results - Struct with fields:
%     sub_g_net     - Selected gene names (K×1 string array)
%     sub_Q_net     - Subnetwork differential co-expression (K×K)
%     sub_Qv        - Per-gene contribution scores (K×1)
%     selected_idx  - Boolean index of selected genes (G×1)
%     G_graph       - Graph object (if plotit=true)
%     Xdiff         - Differential co-expression matrix (G×G)
%     Q1            - Unpenalised QUBO matrix (G×G)
%     Q             - Penalised QUBO matrix (G×G)
%     num_selected  - Number of selected genes
%
% AUTHOR: Selim Romero, Texas A&M University

    % --- Parse options -----------------------------------------------------
    if nargin < 6; opts = struct(); end

    cond_a_label      = getOpt(opts, 'cond_a_label',      'A');
    cond_b_label      = getOpt(opts, 'cond_b_label',      'B');
    method            = getOpt(opts, 'method',            'mnn');
    n_neighbors       = getOpt(opts, 'n_neighbors',       15);
    Vdiff_ext         = getOpt(opts, 'Vdiff',             []);
    cell_state_A      = getOpt(opts, 'cell_state_A',      []);
    cell_state_B      = getOpt(opts, 'cell_state_B',      []);
    use_pathway_prior = getOpt(opts, 'use_pathway_prior', false);
    penalty_scale     = getOpt(opts, 'penalty_scale',     10);
    plotit            = getOpt(opts, 'plotit',            true);
    outfile           = getOpt(opts, 'outfile',           'qubo_genes_sol.txt');
    edge_pct          = getOpt(opts, 'edge_pct',          95);

    batch_id = string(batch_id(:));
    genelist = string(genelist(:));
    g        = string(g(:));

    % ======================================================================
    % Step 1 – Partition X by condition  →  X_A (G × N_A), X_B (G × N_B)
    % ======================================================================
    cs_placeholder = zeros(size(X, 2), 1);
    [X_A, X_B, cs_A_out, cs_B_out] = qubo_dr.preprocess.subset_by_condition( ...
        X, batch_id, cs_placeholder, cond_a_label, cond_b_label);

    fprintf('Condition A (%s): %d cells | Condition B (%s): %d cells\n', ...
        cond_a_label, size(X_A,2), cond_b_label, size(X_B,2));

    % Override cell-state vectors if provided externally
    if ~isempty(cell_state_A); cs_A_out = cell_state_A; end
    if ~isempty(cell_state_B); cs_B_out = cell_state_B; end

    % ======================================================================
    % Step 2 – Non-negative cosine similarity (Gram matrices)
    %          S_A = X_A_norm · X_A_norm^T,  S_B = X_B_norm · X_B_norm^T
    % ======================================================================
    [SA, XA_norm] = qubo_dr.preprocess.compute_gene_similarity(X_A);
    [SB, XB_norm] = qubo_dr.preprocess.compute_gene_similarity(X_B);

    % ======================================================================
    % Step 3 – Differential matrix and V_diff
    %          X_diff = S_B - S_A
    % ======================================================================
    [Xdiff, Vdiff] = qubo_dr.preprocess.compute_differential( ...
        SB, SA, XB_norm, XA_norm, cs_B_out, cs_A_out);

    % Use externally supplied V_diff if provided
    if ~isempty(Vdiff_ext)
        Vdiff = Vdiff_ext(:);
        fprintf('V_diff: using externally supplied vector\n');
    elseif all(cs_A_out == 0) && all(cs_B_out == 0)
        fprintf('V_diff: disabled (no cell-state scalar provided)\n');
    else
        fprintf('V_diff: computed from cell-state scalar projection\n');
    end

    % ======================================================================
    % Step 4 – Subset all matrices to pathway genes
    % ======================================================================
    [~, idx] = ismember(upper(genelist), upper(g));
    idx(idx == 0) = [];
    g_sub       = g(idx);
    Xdiff_sub   = Xdiff(idx, idx);
    Vdiff_sub   = Vdiff(idx);

    % ======================================================================
    % Step 5 – MNN adjacency matrices (gene-space via cell-SVD embedding)
    %          SVD run on cells; gene embedding = U*Sigma in R^(G × d)
    % ======================================================================
    fprintf('Building %s adjacency matrices...\n', upper(method));
    Xmnn_A = qubo_dr.graph.build_adjacency_subset(X_A, method, n_neighbors);
    Xmnn_B = qubo_dr.graph.build_adjacency_subset(X_B, method, n_neighbors);
    Xmnn_A_sub = Xmnn_A(idx, idx);
    Xmnn_B_sub = Xmnn_B(idx, idx);

    % ======================================================================
    % Step 6 – Optional pathway membership prior X_net
    %          Only informative when extra candidate genes are in g_sub.
    % ======================================================================
    Xnet_target = [];
    if use_pathway_prior
        fprintf('Building pathway membership prior X_net...\n');
        Xnet_target = qubo_dr.graph.build_pathway_network(g_sub, genelist);
    else
        fprintf('Pathway prior X_net: skipped (pure pathway-gene analysis)\n');
    end

    % ======================================================================
    % Step 7 – Assemble QUBO matrix
    % ======================================================================
    [Q, Q1, ~] = qubo_dr.qubo.assemble_qubo_matrix( ...
        Xdiff_sub, Vdiff_sub, Xmnn_A_sub, Xmnn_B_sub, K, ...
        Xnet_target, penalty_scale);

    % ======================================================================
    % Step 8 – Solve QUBO
    % ======================================================================
    [selected_idx, sub_g_net, sol] = qubo_dr.qubo.solve_qubo_problem(Q, g_sub);

    % ======================================================================
    % Step 9 – Extract subnetwork
    % ======================================================================
    [sub_Q_net, sub_Qv] = qubo_dr.qubo.extract_subnetwork(Xdiff_sub, Q1, selected_idx, sol);

    % ======================================================================
    % Step 10 – Plots
    % ======================================================================
    G_graph = [];
    if plotit
        figure('Position', [100 100 1200 500]);
        subplot(1,2,1);
        qubo_dr.plot.plot_coexpr_heatmap(sub_Q_net, sub_g_net, ...
            sprintf('Subnetwork Co-expression (K=%d)', K), 'parula');
        subplot(1,2,2);
        [~, G_graph] = qubo_dr.plot.plot_gene_network(sub_Q_net, sub_g_net, edge_pct);
        drawnow;
    end

    % ======================================================================
    % Step 11 – Save results
    % ======================================================================
    write_results(outfile, sub_g_net, sub_Qv);

    % ======================================================================
    % Output struct
    % ======================================================================
    results.sub_g_net    = sub_g_net;
    results.sub_Q_net    = sub_Q_net;
    results.sub_Qv       = sub_Qv;
    results.selected_idx = selected_idx;
    results.G_graph      = G_graph;
    results.Xdiff        = Xdiff_sub;
    results.Q1           = Q1;
    results.Q            = Q;
    results.num_selected = nnz(selected_idx);

    fprintf('Pipeline complete — %d genes selected.\n', results.num_selected);
end

% -------------------------------------------------------------------------
function val = getOpt(opts, field, default)
    if isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default;
    end
end

% -------------------------------------------------------------------------
function write_results(outfile, gene_names, scores)
    fid = fopen(outfile, 'w');
    if fid == -1
        warning('Cannot open %s for writing.', outfile); return;
    end
    fprintf(fid, '%% QUBO differential co-expression — selected subnetwork\n');
    fprintf(fid, '%% Gene\tContribution score\n');
    for i = 1:numel(gene_names)
        fprintf(fid, '%s\t%.6f\n', gene_names(i), scores(i));
    end
    fclose(fid);
end
