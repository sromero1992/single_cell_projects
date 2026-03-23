function scores = compute_ucell_score(X, g, geneset, n_max_rank)
% COMPUTE_UCELL_SCORE  Rank-based gene set activity scoring (UCell method).
%
%   scores = compute_ucell_score(X, g, geneset)
%   scores = compute_ucell_score(X, g, geneset, n_max_rank)
%
%   For each cell, ranks all genes by expression (highest expression =
%   highest rank), then computes the normalized Mann-Whitney U statistic
%   for the genes in 'geneset'. The result is a per-cell activity score
%   in [0, 1] reflecting how highly expressed the gene set is relative to
%   all other genes in that cell.
%
%   INPUT:
%     X          - (G x N) normalized expression matrix (genes x cells)
%     g          - (G x 1) gene name vector (string array or cell array)
%     geneset    - gene names in the pathway/set of interest (string array)
%     n_max_rank - max rank cap (default: 1500). Genes ranked below this
%                  threshold all receive rank n_max_rank, reducing noise
%                  from lowly expressed / dropout genes. Set to G to disable.
%
%   OUTPUT:
%     scores     - (1 x N) per-cell UCell activity score in [0, 1].
%                  Higher values indicate stronger gene set activity.
%
%   ALGORITHM (Andreatta & Carmona, Comput Struct Biotechnol J, 2021):
%     1. For each cell, rank all G genes by expression (rank 1 = lowest,
%        rank G = highest). Ties are broken by average rank.
%     2. Cap ranks at n_max_rank to reduce sensitivity to dropout genes.
%     3. Compute Mann-Whitney U = R - n_s*(n_s+1)/2, where R is the sum
%        of ranks of the gene set genes, and n_s = |geneset|.
%     4. Normalize: score = U / (n_s * (n_max_rank - n_s)).
%
%   REFERENCE:
%     Andreatta M, Carmona SJ. UCell: Robust and scalable single-cell gene
%     signature scoring. Comput Struct Biotechnol J. 2021;19:3796-3798.
%     doi: 10.1016/j.csbj.2021.06.043
%
%   AUTHOR: Selim Romero, Texas A&M University

    % --- Input validation ---
    if nargin < 4 || isempty(n_max_rank)
        n_max_rank = 1500;
    end

    [G, N] = size(X);
    n_max_rank = min(n_max_rank, G);  % Cannot exceed total genes

    % --- Find gene set indices (case-insensitive) ---
    g_upper      = upper(string(g(:)));
    geneset_upper = upper(string(geneset(:)));
    set_mask = ismember(g_upper, geneset_upper);
    n_set    = sum(set_mask);

    if n_set == 0
        warning('qubo_dr:ucell:noGenesFound', ...
            'No genes from geneset found in g. Returning zero scores.');
        scores = zeros(1, N);
        return;
    end

    if n_set >= n_max_rank
        warning('qubo_dr:ucell:setTooLarge', ...
            'Gene set size (%d) >= n_max_rank (%d). Consider increasing n_max_rank.', ...
            n_set, n_max_rank);
    end

    fprintf('  UCell scoring: %d/%d gene set genes found. n_max_rank=%d.\n', ...
        n_set, numel(geneset_upper), n_max_rank);

    % --- Compute ranks and U statistic per cell ---
    scores = zeros(1, N);
    max_U  = n_set * (n_max_rank - n_set);  % Maximum possible U (denominator)

    if max_U <= 0
        warning('qubo_dr:ucell:invalidMaxU', ...
            'max_U <= 0 (n_set=%d, n_max_rank=%d). Returning zeros.', n_set, n_max_rank);
        scores = zeros(1, N);
        return;
    end

    for ci = 1:N
        % Rank genes ascending (rank 1 = lowest expression)
        [~, sort_idx] = sort(X(:, ci), 'ascend');
        ranks = zeros(G, 1);
        ranks(sort_idx) = (1:G)';

        % Cap ranks at n_max_rank (UCell: compress all low-rank genes)
        ranks = min(ranks, n_max_rank);

        % Mann-Whitney U: R = sum of ranks of set genes
        R = sum(ranks(set_mask));

        % U statistic (subtract minimum possible rank sum)
        U = R - n_set * (n_set + 1) / 2;

        % Normalize to [0, 1]
        scores(ci) = U / max_U;
    end

    scores = max(0, min(1, scores));  % Numerical safety clamp
end
