function S = sce_circ_aucell_multi(X, gene_names, gene_sets, top_frac)
%SCE_CIRC_AUCELL_MULTI  AUCell scores for MANY gene sets in one ranking pass.
%   S = sce_circ_aucell_multi(X, gene_names, gene_sets) returns an
%   (nSets x nCell) matrix. Each cell is ranked ONCE and reused for every gene
%   set, so scoring all enriched pathways costs about the same as scoring one.
    if nargin < 4 || isempty(top_frac); top_frac = 0.05; end
    gene_names = upper(string(gene_names));
    [ng, ncell] = size(X);
    np  = numel(gene_sets);
    idx = cell(np,1); nset = zeros(np,1); maxauc = zeros(np,1);
    thr = max(1, ceil(top_frac * ng));
    for p = 1:np
        [tf, loc] = ismember(upper(string(gene_sets{p})), gene_names);
        ii = loc(tf); ii = ii(ii > 0);
        idx{p} = ii; nset(p) = numel(ii);
        maxauc(p) = thr*nset(p) - nset(p)*(nset(p)-1)/2;
        if maxauc(p) <= 0; maxauc(p) = max(1, thr*nset(p)); end
    end
    S = zeros(np, ncell);
    for c = 1:ncell
        x = full(X(:, c));
        [~, ord] = sort(x, 'descend');
        rk = zeros(ng, 1); rk(ord(1:thr)) = 1:thr;        % rank of each top gene
        for p = 1:np
            if nset(p) == 0; continue; end
            hits = rk(idx{p}); hits = hits(hits > 0);
            if ~isempty(hits); S(p, c) = sum(thr - hits + 1) / maxauc(p); end
        end
    end
end
