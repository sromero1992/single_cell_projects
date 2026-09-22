function sce_circ_pathway_record(outdir, celltype, period12, window, term, library, stats, genes)
%SCE_CIRC_PATHWAY_RECORD  Append one pathway's circadian test result to a CSV.
%   Accumulates into  <outdir>/<ct_safe>/<ct_safe>_period_XX_pathway_circadian.csv
%   columns: ZT_window, Term, Library, Acrophase, Amp, Mesor, Pvalue, Corr,
%   Pvalue_corr, Genes.
    if nargin < 3 || isempty(period12); period12 = false; end
    ct_safe = tscp_safe(celltype);
    if period12; per = '_period_12_'; else; per = '_period_24_'; end
    d = fullfile(outdir, ct_safe); if ~exist(d, 'dir'); mkdir(d); end
    f = fullfile(d, [ct_safe per 'pathway_circadian.csv']);
    row = table(string(window), string(term), string(library), stats.acrophase, ...
                stats.amp, stats.mesor, stats.pvalue, stats.corr, stats.pvalue_corr, ...
                string(strjoin(cellstr(genes), ';')), ...
                'VariableNames', {'ZT_window','Term','Library','Acrophase','Amp', ...
                                  'Mesor','Pvalue','Corr','Pvalue_corr','Genes'});
    if isfile(f)
        writetable(row, f, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(row, f);
    end
end
