%% TimeSCape pathway-in-GUI: command-line test (run BEFORE GUI wiring)
% Prereqs in the workspace: sce (your SingleCellExperiment) and tmeta.
% Point the readtable at one cell type's confident circadian list.

ct   = 'CD8+ T cells';                          % cell type label in sce.c_cell_type_tx
csv  = 'CD8_T_cells/CD8_T_cells_period_24_circadian_analysis_confident.csv';
T1c  = readtable(csv);

% 1) bin the circadian list by acrophase (3 h windows; user-definable)
bins = sce_circ_bin_genes(T1c, 3);

% 2) enrich each time bin via Enrichr, keeping a note of the ZT window
lib  = 'KEGG_2019_Mouse';
for i = 1:numel(bins)
    if numel(bins(i).genes) < 5; continue; end
    fprintf('\n== %s : %d genes ==\n', bins(i).label, numel(bins(i).genes));
    bins(i).enrich = sce_circ_enrichr(bins(i).genes, lib);
    if height(bins(i).enrich) > 0
        disp(bins(i).enrich(1:min(5,height(bins(i).enrich)), {'Term','Pvalue','AdjPvalue','nOverlap'}));
    end
end

% 3) pick a bin + pathway, run the circadian test on its AUCell score and plot
sel_bin  = 3;  sel_path = 1;
pw_genes = bins(sel_bin).enrich.Genes{sel_path};
pw_name  = sprintf('%s | %s', bins(sel_bin).label, bins(sel_bin).enrich.Term(sel_path));
figure; ax = axes;
stats = sce_circ_pathway_plot(sce, tmeta, pw_genes, pw_name, ct, false, ax, 'turbo', true);
disp(stats);
