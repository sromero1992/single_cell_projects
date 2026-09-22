function matf = sce_circ_pathway_save(outdir, celltype, period12, bins, library, bin_width)
%SCE_CIRC_PATHWAY_SAVE  Persist binned + enriched pathways for a cell type.
%   Writes into  <outdir>/<ct_safe>/  next to that cell type's circadian CSVs:
%     <ct_safe>_period_XX_pathway_enrichment.mat  (reloadable: bins + config)
%     <ct_safe>_period_XX_pathway_enrichment.csv  (human-readable table)
    if nargin < 3 || isempty(period12); period12 = false; end
    ct_safe = tscp_safe(celltype);
    if period12; per = '_period_12_'; else; per = '_period_24_'; end
    d = fullfile(outdir, ct_safe);
    if ~exist(d, 'dir'); mkdir(d); end
    info = struct('celltype', char(celltype), 'library', char(library), ...
                  'bin_width', bin_width, 'period12', period12, ...
                  'saved', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'bins', bins);
    matf = fullfile(d, [ct_safe per 'pathway_enrichment.mat']);
    save(matf, 'info');
    rows = {};
    for i = 1:numel(bins)
        if ~isfield(bins, 'enrich') || isempty(bins(i).enrich); continue; end
        E = bins(i).enrich;
        for k = 1:height(E)
            rows(end+1, :) = { bins(i).label, char(E.Term(k)), E.Pvalue(k), ...
                E.AdjPvalue(k), E.nOverlap(k), strjoin(E.Genes{k}, ';') }; %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        T = cell2table(rows, 'VariableNames', ...
              {'ZT_window','Term','Pvalue','AdjPvalue','nOverlap','Genes'});
        writetable(T, fullfile(d, [ct_safe per 'pathway_enrichment.csv']));
    end
end
