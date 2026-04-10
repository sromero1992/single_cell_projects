function results = bootstrap_stability(X, g, batch_id, genelist, K, opts)
% BOOTSTRAP_STABILITY  Bootstrap cell resampling to assess gene selection
%                      stability in the QUBO subnetwork.
%
%   results = bootstrap_stability(X, g, batch_id, genelist, K, opts)
%
%   For each bootstrap iteration, cells are resampled with replacement
%   within each condition. The full QUBO pipeline (similarity matrices,
%   MNN adjacency, QUBO assembly, QUBO solve) is re-run on the resampled
%   data. Gene selection frequencies across bootstraps quantify which
%   genes are robustly selected versus incidentally selected.
%
%   INPUT:
%     X        - (G x N) normalized expression matrix (pathway genes only)
%     g        - (G x 1) gene names (string array)
%     batch_id - (1 x N) condition label per cell (string array)
%     genelist - pathway gene list (string array), used to build Xnet_target
%     K        - target subnetwork size
%     opts     - struct with optional fields:
%       .condition_A  - label for condition A (default: 'KO')
%       .condition_B  - label for condition B (default: 'WT')
%       .n_boot       - number of bootstrap iterations (default: 200)
%       .method       - neighbor method: 'mnn' or 'knn' (default: 'mnn')
%       .n_neighbors  - K for MNN graph (default: 30)
%       .seed         - random seed (default: 42)
%       .verbose      - print progress (default: true)
%       .freq_thresh  - frequency threshold for "stable" genes (default: 0.5)
%
%   OUTPUT:
%     results - struct with fields:
%       .gene_frequency  - (G x 1) selection frequency per gene in [0,1]
%       .stability_index - scalar: mean frequency of genes above freq_thresh
%       .freq_matrix     - (G x n_boot) binary selection matrix
%       .stable_genes    - gene names selected in > freq_thresh of boots
%       .gene_names      - gene names in order (same order as gene_frequency)
%       .freq_thresh     - threshold used
%
%   INTERPRETATION:
%     genes with gene_frequency > 0.5 are "core stable" — they are selected
%     in the majority of bootstrap samples, indicating the QUBO solution is
%     robust to cell sampling variability. stability_index near 1.0 means
%     the full subnetwork is highly reproducible.
%
%   AUTHOR: Selim Romero, Texas A&M University

    % --- Default options ---
    if nargin < 6 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'condition_A');  opts.condition_A = 'KO';   end
    if ~isfield(opts, 'condition_B');  opts.condition_B = 'WT';   end
    if ~isfield(opts, 'n_boot');       opts.n_boot = 200;         end
    if ~isfield(opts, 'method');       opts.method = 'mnn';       end
    if ~isfield(opts, 'n_neighbors'); opts.n_neighbors = 30;     end
    if ~isfield(opts, 'seed');         opts.seed = 42;            end
    if ~isfield(opts, 'verbose');      opts.verbose = true;       end
    if ~isfield(opts, 'freq_thresh'); opts.freq_thresh = 0.5;    end

    rng(opts.seed);

    [G, N] = size(X);

    % --- Condition indices ---
    idx_A = find(contains(batch_id, opts.condition_A));
    idx_B = find(contains(batch_id, opts.condition_B));

    if isempty(idx_A) || isempty(idx_B)
        error('qubo_dr:bootstrap:noConditionCells', ...
            'No cells found for condition A ("%s") or B ("%s").', ...
            opts.condition_A, opts.condition_B);
    end

    % --- Pathway prior matrix (static across bootstraps) ---
    idx_pathway = ismember(upper(g), upper(string(genelist(:))));
    Xnet = zeros(G, G);
    for ig = 1:G
        if idx_pathway(ig)
            Xnet(ig, idx_pathway) = -1;
            Xnet(idx_pathway, ig) = -1;
        end
    end

    % --- Bootstrap loop ---
    freq_matrix = zeros(G, opts.n_boot);

    if opts.verbose
        fprintf('\nRunning bootstrap stability (%d iterations)...\n', opts.n_boot);
    end

    n_successful = 0;
    for b = 1:opts.n_boot

        if opts.verbose && mod(b, 50) == 0
            fprintf('  Bootstrap %d / %d  (successful: %d)\n', b, opts.n_boot, n_successful);
        end

        try
            % --- Resample cells with replacement within each condition ---
            boot_A_idx = idx_A(randi(numel(idx_A), numel(idx_A), 1));
            boot_B_idx = idx_B(randi(numel(idx_B), numel(idx_B), 1));

            Xa_b = X(:, boot_A_idx);
            Xb_b = X(:, boot_B_idx);

            % --- Similarity matrices ---
            Za_b   = normr(Xa_b);
            Zb_b   = normr(Xb_b);
            Sa_b   = Za_b * Za_b';
            Sb_b   = Zb_b * Zb_b';
            dS_b = Sb_b - Sa_b;

            % --- MNN adjacency matrices ---
            Amnn_A = adjX_mat_construct_sparse_cust_idx(Xa_b, opts.method, opts.n_neighbors);
            Amnn_B = adjX_mat_construct_sparse_cust_idx(Xb_b, opts.method, opts.n_neighbors);

            % --- QUBO assembly ---
            Q1_b = dS_b - full(Amnn_A + Amnn_B) + Xnet;
            P    = 10 * max(abs(Q1_b(:)));

            Q_b = Q1_b + P * (ones(G, G) - eye(G));
            for i = 1:G
                Q_b(i, i) = Q_b(i, i) + P * (1 - 2 * K);
            end

            % --- Solve QUBO ---
            qp  = qubo(Q_b);
            res = solve(qp);

            freq_matrix(:, b) = double(res.BestX == 1);
            n_successful = n_successful + 1;

        catch ME
            if opts.verbose
                fprintf('  [Warning] Bootstrap %d failed: %s\n', b, ME.message);
            end
            % Leave freq_matrix(:, b) = 0 (gene not selected in failed boot)
        end
    end

    if n_successful == 0
        error('qubo_dr:bootstrap:allFailed', ...
            'All %d bootstrap iterations failed. Check inputs.', opts.n_boot);
    end

    % --- Compute frequencies (over successful boots only) ---
    gene_frequency = sum(freq_matrix, 2) / n_successful;  % (G x 1)

    % --- Stable genes ---
    stable_mask   = gene_frequency > opts.freq_thresh;
    stable_genes  = g(stable_mask);

    if any(stable_mask)
        stability_index = mean(gene_frequency(stable_mask));
    else
        stability_index = 0.0;
    end

    % --- Pack results ---
    results.gene_frequency  = gene_frequency;
    results.stability_index = stability_index;
    results.freq_matrix     = freq_matrix;
    results.stable_genes    = stable_genes;
    results.gene_names      = g;
    results.freq_thresh     = opts.freq_thresh;
    results.n_successful    = n_successful;

    % --- Report ---
    if opts.verbose
        fprintf('\n--- Bootstrap Stability Results (%d/%d successful) ---\n', ...
            n_successful, opts.n_boot);
        fprintf('  Stable genes (freq > %.2f): %d\n', opts.freq_thresh, sum(stable_mask));
        fprintf('  Stability index:            %.3f\n', stability_index);
        if ~isempty(stable_genes)
            top_n = min(10, numel(stable_genes));
            fprintf('  Top stable genes:           %s\n', ...
                strjoin(stable_genes(1:top_n), ', '));
        end
        fprintf('------------------------------------------------------\n');
    end
end
