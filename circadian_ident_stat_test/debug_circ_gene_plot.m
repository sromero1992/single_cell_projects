% Debugging fitting
% Celltype
ic = find(sce.c_cell_type_tx == 'T cells');
sce_sub = sce.selectcells(ic);
sce_sub = sce_sub.qcfilter;
sce_sub.X = sc_norm(full(sce_sub.X));

ig = find(sce_sub.g == 'Capn5');

batch_time = unique(sce_sub.c_batch_id);
nzts = length(batch_time);

Xg_zts = {};
% Prepare gene expression per time point
for it = 1:nzts
    ics = find(sce_sub.c_batch_id == batch_time(it));
    % Gene expression for a time point
    if ~isempty(ig)
        Xg_zts{it} = full(sce_sub.X(ig, ics));
    else
        Xg_zts{it} = zeros(size(ics));
    end
end

[tmp_acro, tmp_amp, tmp_T, tmp_mesor, tmp_p_value] = ...
    estimate_phaseR(Xg_zts, 3, false, 'Ftest')


for i = 1:8 
    mean(Xg_zts{i})
end