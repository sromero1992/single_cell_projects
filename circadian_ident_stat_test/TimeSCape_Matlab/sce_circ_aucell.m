function scores = sce_circ_aucell(X, gene_names, gene_set, top_frac)
%SCE_CIRC_AUCELL  AUCell-style pathway activity score per cell.
%   Mirrors Bioconductor AUCell: for each cell, rank genes by expression and
%   compute the area under the recovery curve of GENE_SET within the top
%   TOP_FRAC of the ranking. Depth-robust; returns a 1 x Ncell vector in [0,1].
%
%   X          - genes x cells matrix (sparse or full), raw or normalized.
%   gene_names - string/cellstr of length size(X,1).
%   gene_set   - string/cellstr of pathway gene symbols.
%   top_frac   - ranking threshold (default 0.05 = top 5%).
    if nargin < 4 || isempty(top_frac); top_frac = 0.05; end
    gene_names = string(gene_names);
    gene_set   = string(gene_set);
    [ng, ncell] = size(X);
    inset = ismember(upper(gene_names), upper(gene_set));
    nset  = sum(inset);
    scores = zeros(1, ncell);
    if nset == 0
        warning('sce_circ_aucell: none of the pathway genes are in the data.');
        return;
    end
    thr = max(1, ceil(top_frac * ng));
    maxauc = thr * nset - nset * (nset - 1) / 2;   % area if all set genes ranked top
    if maxauc <= 0; maxauc = thr * nset; end
    for c = 1:ncell
        x = full(X(:, c));
        [~, ord] = sort(x, 'descend');
        hits = find(inset(ord(1:thr)));          % ranks where set genes appear
        if ~isempty(hits)
            scores(c) = sum(thr - hits + 1) / maxauc;   % normalized AUC in [0,1]
        end
    end
end
