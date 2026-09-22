function T = sce_circ_pathway_all_rhythms(sce, tmeta, bins, period12, cust_cells, norm_str, outdir)
%SCE_CIRC_PATHWAY_ALL_RHYTHMS  Circadian test for EVERY enriched pathway.
%   Scores all enriched pathways (shared ranking) and cosinor-tests each,
%   writing <ct>/<ct>_period_XX_pathway_circadian_all.csv with the enrichment
%   stats AND the rhythm stats (Acrophase, Amp, Mesor, p, Corr, corr p).
    if nargin < 4 || isempty(period12); period12 = false; end
    if nargin < 6 || isempty(norm_str); norm_str = 'lib_size'; end
    if nargin < 7 || isempty(outdir);   outdir   = pwd; end
    if period12; period = 12; else; period = 24; end

    ic0 = find(sce.c_cell_type_tx == cust_cells);
    if isempty(ic0); error('No cells of type "%s".', cust_cells); end
    Xsub = sce.X(:, ic0); batch_sub = sce.c_batch_id(ic0);
    if strcmp(norm_str, 'lib_size'); Xn = log1p(pkg.norm_libsize(Xsub, 1e4));
    else;                            Xn = Xsub; end

    rows = {}; gene_sets = {};
    for i = 1:numel(bins)
        if ~isfield(bins,'enrich') || isempty(bins(i).enrich); continue; end
        E = bins(i).enrich;
        for k = 1:height(E)
            gene_sets{end+1} = E.Genes{k}; %#ok<AGROW>
            rows(end+1,:) = { bins(i).label, char(E.Term(k)), E.Pvalue(k), ...
                              E.AdjPvalue(k), E.nOverlap(k), strjoin(E.Genes{k}, ';') }; %#ok<AGROW>
        end
    end
    if isempty(gene_sets); T = table(); return; end

    S = sce_circ_aucell_multi(Xn, sce.g, gene_sets);       % nPathway x nCell

    batch_time = unique(batch_sub); nz = numel(batch_time); at = nan(nz,1);
    for it = 1:nz
        kk = find(tmeta.new_labels == batch_time(it), 1);
        if ~isempty(kk); at(it) = tmeta.ZT_times(kk); end
    end
    valid = ~isnan(at); batch_time = batch_time(valid); at = at(valid);
    [at, ord] = sort(at); batch_time = batch_time(ord); nz = numel(at);

    np = size(S,1);
    Acro=zeros(np,1); Amp=zeros(np,1); Mesor=zeros(np,1); Pv=zeros(np,1); Corr=zeros(np,1); Pc=zeros(np,1);
    for p = 1:np
        Xg = cell(1, nz);
        for it = 1:nz; ics = find(batch_sub == batch_time(it)); Xg{it} = S(p, ics); end
        [ac, am, ~, me, pv, rh, pc] = estimate_phaseR(Xg, at', period12, 'Ftest');
        Acro(p)=ac; Amp(p)=am; Mesor(p)=me; Pv(p)=pv; Corr(p)=rh; Pc(p)=pc;
    end
    meta = cell2table(rows, 'VariableNames', ...
        {'ZT_window','Term','Enrich_Pvalue','Enrich_AdjPvalue','nOverlap','Genes'});
    T = [ meta(:, {'ZT_window','Term','Enrich_Pvalue','Enrich_AdjPvalue','nOverlap'}), ...
          table(Acro, Amp, Mesor, Pv, Corr, Pc, 'VariableNames', ...
                {'Acrophase','Amp','Mesor','Rhythm_Pvalue','Corr','Rhythm_Pvalue_corr'}), ...
          meta(:, 'Genes') ];
    T = sortrows(T, 'Rhythm_Pvalue', 'ascend');

    ct_safe = tscp_safe(cust_cells);
    if period12; per = '_period_12_'; else; per = '_period_24_'; end
    d = fullfile(outdir, ct_safe); if ~exist(d,'dir'); mkdir(d); end
    writetable(T, fullfile(d, [ct_safe per 'pathway_circadian_all.csv']));
end
