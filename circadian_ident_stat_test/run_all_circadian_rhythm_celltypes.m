times = [0 3 6 9 12 15 18 21]';
old_labels = unique(sce.c_batch_id);
% Convert times to ZT labels
new_labels = cell(numel(times), 1);
for i = 1:numel(times)
    if times(i) >= 0
        new_labels{i} = sprintf('ZT%02d', times(i));
    else
        new_labels{i} = sprintf('%d', times(i));
    end
end

new_labels = string(new_labels);
tmeta = table(old_labels, new_labels, times);
[T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta);