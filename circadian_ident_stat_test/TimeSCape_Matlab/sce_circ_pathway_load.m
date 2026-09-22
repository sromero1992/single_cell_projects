function info = sce_circ_pathway_load(outdir, celltype, period12)
%SCE_CIRC_PATHWAY_LOAD  Load saved pathway enrichment for a cell type ([] if none).
%   Returns INFO with fields bins, library, bin_width, celltype, saved — so you
%   can re-plot pathways as post-processing without re-running Enrichr.
    if nargin < 3 || isempty(period12); period12 = false; end
    ct_safe = tscp_safe(celltype);
    if period12; per = '_period_12_'; else; per = '_period_24_'; end
    matf = fullfile(outdir, ct_safe, [ct_safe per 'pathway_enrichment.mat']);
    if ~isfile(matf); info = []; return; end
    S = load(matf, 'info'); info = S.info;
end
