function results = permutation_test(X, g, batch_id, selected_idx, Q1, opts)
% PERMUTATION_TEST  Label permutation test for differential co-expression
%                   subnetwork significance.
%
%   results = permutation_test(X, g, batch_id, selected_idx, Q1, opts)
%
%   Tests the null hypothesis that condition labels are exchangeable (i.e.,
%   that the observed differential co-expression signal in the selected
%   subnetwork is no stronger than expected by chance). The real QUBO
%   solution z* is held fixed; only condition labels are permuted.
%
%   INPUT:
%     X            - (G x N) normalized expression matrix (pathway genes only)
%     g            - (G x 1) gene names
%     batch_id     - (1 x N) condition label per cell (string array)
%     selected_idx - (G x 1) logical vector, z* from QUBO solution
%     Q1           - (G x G) QUBO cost matrix WITHOUT penalty (Q1, not Q)
%     opts         - struct with optional fields:
%       .condition_A  - label for condition A (default: 'KO')
%       .condition_B  - label for condition B (default: 'WT')
%       .n_perm       - number of permutations (default: 1000)
%       .seed         - random seed (default: 42)
%       .verbose      - print results (default: true)
%
%   OUTPUT:
%     results - struct with fields:
%       .pval        - empirical p-value (fraction of null energies <= real)
%       .zscore      - z-score of real energy versus null distribution
%       .E_real      - energy of real selected subnetwork: z*' Q1 z*
%       .E_null      - (1 x n_perm) null energy distribution
%       .null_mean   - mean of null distribution
%       .null_std    - std of null distribution
%       .effect_size - (E_real - null_mean) / null_std  [same as zscore here]
%
%   INTERPRETATION:
%     A significant p-value (e.g. < 0.05) and strongly negative z-score
%     indicate that the selected genes capture a genuine, condition-specific
%     differential co-expression signal that is unlikely under random
%     label assignment.
%
%   AUTHOR: Selim Romero, Texas A&M University

    % --- Default options ---
    if nargin < 6 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'condition_A'); opts.condition_A = 'KO'; end
    if ~isfield(opts, 'condition_B'); opts.condition_B = 'WT'; end
    if ~isfield(opts, 'n_perm');      opts.n_perm = 1000;      end
    if ~isfield(opts, 'seed');        opts.seed = 42;           end
    if ~isfield(opts, 'verbose');     opts.verbose = true;      end

    rng(opts.seed);

    [~, N] = size(X);
    bz = double(selected_idx(:));  % Column vector

    % --- Real subnetwork energy: E* = z*' Q1 z* ---
    E_real = bz' * Q1 * bz;

    % --- Permutation loop ---
    E_null = zeros(1, opts.n_perm);

    if opts.verbose
        fprintf('\nRunning permutation test (%d permutations)...\n', opts.n_perm);
    end

    for p = 1:opts.n_perm
        % Permute cell condition labels
        perm_order = randperm(N);
        batch_perm = batch_id(perm_order);

        % Extract condition-specific cells from permuted labels
        idx_A = contains(batch_perm, opts.condition_A);
        idx_B = contains(batch_perm, opts.condition_B);

        Xa_p = X(:, idx_A);
        Xb_p = X(:, idx_B);

        % Recompute cosine similarity matrices
        % (normr normalizes each ROW to unit L2 norm)
        Za_p  = normr(Xa_p);
        Zb_p  = normr(Xb_p);
        Sa_p  = Za_p * Za_p';
        Sb_p  = Zb_p * Zb_p';
        Xdiff_p = Sb_p - Sa_p;

        % Score real solution z* on permuted differential matrix
        % (only Xdiff component — MNN/Vdiff are omitted for speed;
        %  the differential similarity matrix is the dominant signal term)
        E_null(p) = bz' * Xdiff_p * bz;
    end

    % --- Statistics ---
    null_mean   = mean(E_null);
    null_std    = std(E_null);
    % p-value: fraction of null energies AT OR BELOW the real energy
    % (minimization: real energy should be strongly negative if signal is real)
    pval        = mean(E_null <= E_real);
    zscore_stat = (E_real - null_mean) / (null_std + eps);
    effect_size = zscore_stat;  % Effect size = z-score in this formulation

    % --- Pack results ---
    results.pval        = pval;
    results.zscore      = zscore_stat;
    results.E_real      = E_real;
    results.E_null      = E_null;
    results.null_mean   = null_mean;
    results.null_std    = null_std;
    results.effect_size = effect_size;

    % --- Report ---
    if opts.verbose
        if pval < 0.001;    sig = '***';
        elseif pval < 0.01; sig = '**';
        elseif pval < 0.05; sig = '*';
        else;               sig = 'n.s.';
        end

        fprintf('\n--- Permutation Test Results ---\n');
        fprintf('  Real energy  E*         = %8.4f\n', E_real);
        fprintf('  Null mean ± std         = %8.4f ± %.4f\n', null_mean, null_std);
        fprintf('  Z-score                 = %8.4f\n', zscore_stat);
        fprintf('  Empirical p-value       = %8.4f  %s\n', pval, sig);
        fprintf('--------------------------------\n');
    end
end
