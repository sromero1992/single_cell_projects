function [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, rm_low_conf, period12, ...
                                    custom_genelist, custom_celltype, plot_heat, norm_str, num_cores)
% SCE_CIRC_PHASE_ESTIMATION_STATTEST
%   TimeSCape circadian-rhythm detection pipeline on a SingleCellExperiment.
%
% KEY ROBUSTNESS FEATURES:
%   • Missing time points: if a cell type lacks cells at one or more ZT
%     times, the cosine is fitted using only the time points that ARE
%     present (at their true ZT hours).  No imputation, no zeros.
%   • Parallel block processing for large gene sets.
%   • Dual output: ALL genes + a second confident-only file.
%
% USAGE:
%   times      = [0 3 6 9 12 15 18 21]';
%   old_labels = unique(sce.c_batch_id);
%   new_labels = string(arrayfun(@(t) sprintf('ZT%02d',t), times, 'UniformOutput',false));
%   tmeta      = table(old_labels, new_labels, times, 'VariableNames',{'old_labels','new_labels','ZT_times'});
%   [T1, T2]   = sce_circ_phase_estimation_stattest(sce, tmeta);
%
% INPUTS:
%   sce            - SingleCellExperiment object.
%   tmeta          - Table: old_labels | new_labels | ZT_times (numeric hrs).
%   rm_low_conf    - (default true) Write the confident-only secondary files.
%   period12       - (default false) Use 12-hr period.
%   custom_genelist- (default []) Restrict to these genes.
%   custom_celltype- (default []) Restrict to this cell type.
%   plot_heat      - (default true) Generate heatmap after analysis.
%   norm_str       - (default 'lib_size') 'lib_size' | 'none' | 'magic_impute'.
%   num_cores      - (default max(2,ceil(numcores/4))) parallel pool workers.
%
% OUTPUTS:
%   T1  - Circadian statistics table, ALL genes (no p-value filter).
%   T2  - Per-timepoint mean expression (rows match T1).
%
% SAVED FILES (per cell type, prefix = "<CellType>_<period>_"):
%   *circadian_analysis_all.csv              - All genes (T1)
%   *circadian_analysis_confident.csv        - Both-test p<0.05 genes
%   *circadian_ZTs_mean.csv                  - Raw ZT means (T2, all)
%   *circadian_ZTs_mean_normalized.csv       - ZT00-normalised (T3, all)
%   *circadian_ZTs_mean_confident.csv        - ZT means, confident only
%   *circadian_ZTs_mean_normalized_confident.csv
%   *summary_results.csv                     - Cross-cell-type summary

    tic;
    rng('default');

    % ── Defaults ───────────────────────────────────────────────────────────
    if nargin < 3  || isempty(rm_low_conf);     rm_low_conf     = true;       end
    if nargin < 4  || isempty(period12);         period12        = false;      end
    if nargin < 5  || isempty(custom_genelist);  custom_genelist = {};         end
    if nargin < 6  || isempty(custom_celltype);  custom_celltype = {};         end
    if nargin < 7  || isempty(plot_heat);        plot_heat       = true;       end
    if nargin < 8  || isempty(norm_str);         norm_str        = 'lib_size'; end
    if nargin < 9  || isempty(num_cores)
        num_cores = max(2, ceil(feature('numcores') / 4));
    end
    num_cores = max(1, round(num_cores));

    if period12; per_label = "_period_12_"; else; per_label = "_period_24_"; end
    disp("Circadian identification — period: " + strtrim(per_label));

    % ── Re-label batches to ZT names ───────────────────────────────────────
    batches    = unique(sce.c_batch_id);
    %tmeta.times = sortrows(tmeta.times);
    tmeta.ZT_times = sortrows(tmeta.ZT_times);

    for ib = 1:length(batches)
        str_idx = find(batches(ib) == tmeta.old_labels);
        if isempty(str_idx); continue; end
        idx = find(sce.c_batch_id == batches(ib));
        sce.c_batch_id(idx) = tmeta.new_labels(str_idx);
    end
    batches = unique(sce.c_batch_id);
    disp("  Batches after re-labelling: " + strjoin(batches', ', '));

    % ── Parallel pool ──────────────────────────────────────────────────────
    fprintf('  Parallel workers requested: %d\n', num_cores);
    pool = gcp('nocreate');
    if isempty(pool)
        % No pool exists — start fresh
        parpool(num_cores);
        tic_warm = tic;
        parfor i = 1:num_cores; pause(0.01); end  % warm up workers
        fprintf('  Parallel pool started and warmed up (%d workers, %.1f s)\n', ...
                num_cores, toc(tic_warm));
    elseif pool.NumWorkers ~= num_cores
        % Pool exists but wrong size — restart with correct count
        fprintf('  Resizing pool: %d → %d workers…\n', pool.NumWorkers, num_cores);
        delete(pool);
        parpool(num_cores);
        tic_warm = tic;
        parfor i = 1:num_cores; pause(0.01); end
        fprintf('  Parallel pool resized and warmed up (%d workers, %.1f s)\n', ...
                num_cores, toc(tic_warm));
    else
        % Correct pool already running — reuse, no warm-up needed
        fprintf('  Parallel pool already running (%d workers) — reusing.\n', pool.NumWorkers);
    end

    % ── Normalisation strategy ─────────────────────────────────────────────
    % MEMORY NOTE: For large datasets (>500k cells) 'lib_size' normalises
    % per-cell-type inside the loop below to avoid holding the full dense
    % matrix in RAM.  For lib-size/10k this is mathematically IDENTICAL to
    % normalising the full matrix first, because each cell's size factor
    % depends only on its own total count sum — not on which other cells
    % or genes are present.  Cancer-stage mixing, replicate mixing, etc.
    % do NOT affect the result.
    %
    % norm_str options:
    %   'lib_size'     — library-size to 10k + log1p  (default, recommended)
    %   'none'         — use sce.X as-is (pass already-normalised data,
    %                    or raw counts if you know what you are doing)
    %   'magic_impute' — MAGIC imputation (slow, for cancer/dropout-heavy data)
    fprintf('  Normalisation: %s\n', norm_str);

    % For MAGIC we must normalise the full object upfront (MAGIC is not
    % per-cell-independent).  Store in X_magic — never mutate sce.X so
    % that repeated per-cell-type calls all see the original handle object.
    X_magic = [];
    if strcmp(norm_str, 'magic_impute')
        X_magic = sparse(sc_impute(sce.X, 'MAGIC'));
    end
    % lib_size and none: sce.X left untouched; normalisation (or not)
    % happens inside the cell-type loop to keep peak RAM minimal.

    % ── Cell-type loop ─────────────────────────────────────────────────────
    cell_type_list = unique(sce.c_cell_type_tx);
    ncell_types    = length(cell_type_list);
    nztps_expected = length(unique(tmeta.new_labels));  % expected # of ZT labels

    info_p_type = zeros(ncell_types, 9);

    for icell_type = 1:ncell_types

        cell_type = cell_type_list(icell_type);
        if ~isempty(custom_celltype)
            if ~ismember(cell_type, custom_celltype); continue; end
        end
        fprintf("\nProcessing cell type: %s\n", cell_type);

        % ── Create per-cell-type output directory ──────────────────────────
        outdir_name = regexprep(strtrim(char(cell_type)), '[^a-zA-Z0-9_]', '_');
        outdir_name = regexprep(outdir_name, '_+', '_');   % collapse runs: Cd163__Mrc1 → Cd163_Mrc1
        outdir_name = regexprep(outdir_name, '^_|_$', ''); % strip leading/trailing underscores
        outdir      = fullfile(pwd, outdir_name);
        if ~exist(outdir, 'dir'); mkdir(outdir); end
        fprintf("  Output directory: %s\n", outdir);

        % ── Direct-index subset — NEVER call selectcells() on a handle object ─
        % selectcells() mutates the caller's sce in place (handle semantics),
        % permanently restricting sce.c_cell_type_tx for all future iterations.
        % Instead, keep all data in sce untouched and use column indices.
        ct_idx     = find(sce.c_cell_type_tx == cell_type);
        n_cells_ct = length(ct_idx);

        % ── Per-cell-type normalisation ────────────────────────────────────
        % lib_size: peak RAM = one dense submatrix at a time; result is
        % bit-for-bit identical to normalising the full matrix because
        % lib-size depends only on each cell's own total count sum.
        if strcmp(norm_str, 'lib_size')
            X_sub = pkg.norm_libsize(sce.X(:, ct_idx), 1e4);
            X_sub = sparse(log1p(X_sub));
        elseif strcmp(norm_str, 'magic_impute')
            X_sub = X_magic(:, ct_idx);
        else  % 'none'
            X_sub = sce.X(:, ct_idx);
        end
        batch_sub = sce.c_batch_id(ct_idx);
        clear ct_idx;

        % Gene list
        if isempty(custom_genelist)
            gene_list = sce.g;
        else
            if ischar(custom_genelist) || isstring(custom_genelist)
                custom_genelist = cellstr(custom_genelist);
            end
            if iscolumn(custom_genelist); custom_genelist = custom_genelist'; end
            [lic, ~] = ismember(custom_genelist, sce.g);
            gene_list = custom_genelist(lic);
            disp("  Matching genes: " + strjoin(gene_list, ', '));
        end
        num_genes = length(gene_list);

        % ── Time-point discovery ───────────────────────────────────────────
        % Only use time points that actually have cells.
        % Build actual_times: the true ZT hour for each present batch.
        batch_time = unique(batch_sub);
        nzts       = length(batch_time);

        actual_times = nan(nzts, 1);
        for it = 1:nzts
            idx_t = find(tmeta.new_labels == batch_time(it), 1);
            if ~isempty(idx_t)
                %actual_times(it) = tmeta.times(idx_t);
                actual_times(it) = tmeta.ZT_times(idx_t);
            end
        end

        % Drop any time points we could not map to a numeric ZT hour
        valid_tp  = ~isnan(actual_times);
        batch_time   = batch_time(valid_tp);
        actual_times = actual_times(valid_tp);
        nzts         = sum(valid_tp);

        fprintf("  Time points found: %d / %d expected\n", nzts, nztps_expected);
        if nzts < 4
            warning("  Cell type '%s' has fewer than 4 time points — skipping.", cell_type);
            continue;
        end
        if nzts < nztps_expected
            fprintf("  ⚠ Missing %d time point(s) — fitting on available points only.\n", ...
                    nztps_expected - nzts);
        end

        % ── Block-parallel cosinor fitting ─────────────────────────────────
        if n_cells_ct > 15000; block_size = 50; else; block_size = 500; end
        num_blocks = ceil(num_genes / block_size);

        acro        = zeros(num_genes, 1);
        amp         = zeros(num_genes, 1);
        T           = zeros(num_genes, 1);
        mesor_arr   = zeros(num_genes, 1);
        R0          = zeros(num_genes, nzts);
        p_value     = zeros(num_genes, 1);
        rho         = zeros(num_genes, 1);
        p_value_rho = zeros(num_genes, 1);

        % Snapshot variables needed inside parfor (plain arrays, not handles)
        actual_times_snap = actual_times;
        batch_time_snap   = batch_time;
        sce_sub_X         = X_sub;
        sce_sub_batch     = batch_sub;
        sce_sub_g         = sce.g;

        progressbar('Starting circadian analysis...');

        for block_idx = 1:num_blocks
            s_idx = (block_idx - 1) * block_size + 1;
            e_idx = min(block_idx * block_size, num_genes);
            cur_block = gene_list(s_idx:e_idx);
            n_blk     = length(cur_block);

            tmp_acro        = zeros(n_blk, 1);
            tmp_amp         = zeros(n_blk, 1);
            tmp_T           = zeros(n_blk, 1);
            tmp_p_value     = zeros(n_blk, 1);
            tmp_R0          = zeros(n_blk, nzts);
            tmp_mesor       = zeros(n_blk, 1);
            tmp_rho         = zeros(n_blk, 1);
            tmp_p_value_rho = zeros(n_blk, 1);

            gene_indices = zeros(n_blk, 1);
            for gi = 1:n_blk
                gene_indices(gi) = find(sce_sub_g == string(cur_block{gi}), 1);
            end

            parfor igene_blk = 1:n_blk
                ig     = gene_indices(igene_blk);
                Xg_zts = cell(1, nzts);

                for it = 1:nzts
                    ics = find(sce_sub_batch == batch_time_snap(it));
                    if ~isempty(ig) && ~isempty(ics)
                        Xg_zts{it}              = full(sce_sub_X(ig, ics));
                        tmp_R0(igene_blk, it)   = mean(Xg_zts{it});
                    else
                        Xg_zts{it}              = [];   % empty → skipped inside estimate_phaseR
                        tmp_R0(igene_blk, it)   = NaN;
                    end
                end

                [tmp_acro(igene_blk), tmp_amp(igene_blk), tmp_T(igene_blk), ...
                 tmp_mesor(igene_blk), tmp_p_value(igene_blk), ...
                 tmp_rho(igene_blk), tmp_p_value_rho(igene_blk)] = ...
                    estimate_phaseR(Xg_zts, actual_times_snap', period12, 'Ftest');
            end

            acro(s_idx:e_idx)        = tmp_acro;
            amp(s_idx:e_idx)         = tmp_amp;
            T(s_idx:e_idx)           = tmp_T;
            mesor_arr(s_idx:e_idx)   = tmp_mesor;
            R0(s_idx:e_idx, :)       = tmp_R0;
            p_value(s_idx:e_idx)     = tmp_p_value;
            rho(s_idx:e_idx)         = tmp_rho;
            p_value_rho(s_idx:e_idx) = tmp_p_value_rho;

            clear tmp_acro tmp_amp tmp_T tmp_mesor tmp_R0 ...
                  tmp_p_value tmp_rho tmp_p_value_rho cur_block gene_indices;

            progressbar((block_idx / num_blocks) * 100);
        end

        % ── Wrap acrophase into [0, 24) ────────────────────────────────────
        acro_fmt            = acro;
        acro_fmt(acro <  0) = acro(acro <  0) + 24;
        acro_fmt(acro > 24) = acro(acro > 24) - 24;

        % ── BH-adjusted p-values ───────────────────────────────────────────
        p_adj     = bh_adjust_pvalues(p_value);
        p_adj_rho = bh_adjust_pvalues(p_value_rho);

        % ── Result tables ──────────────────────────────────────────────────
        T1 = table(gene_list, amp, abs(amp), mesor_arr, acro, acro_fmt, T, ...
                   p_value, p_adj, rho, p_value_rho, p_adj_rho);
        T1.Properties.VariableNames = ["Genes","Amp","Abs_Amp","Mesor", ...
                                        "Acrophase","Acrophase_24","Period", ...
                                        "pvalue","pvalue_adj", ...
                                        "Sine_corr","pvalue_corr","pvalue_adj_corr"];

        % T2: dynamic columns from actual time labels
        zt_labels   = string(batch_time(:)');
        T2_varnames = ["Genes", zt_labels];
        %T2          = cell2table([gene_list, num2cell(R0)], 'VariableNames', T2_varnames);
        T2 = array2table(R0, 'VariableNames', zt_labels);
        T2 = [table(gene_list, 'VariableNames', {'Genes'}), T2];

        % Remove genes where either p-value is NaN (fit failed)
        valid = ~isnan(T1.pvalue) & ~isnan(T1.pvalue_corr);
        T1 = T1(valid, :);
        T2 = T2(valid, :);

        % ── Confidence statistics ──────────────────────────────────────────
        conf_ftest = T1.pvalue      < 0.05;
        conf_corr  = T1.pvalue_corr < 0.05;
        conf_both  = conf_ftest & conf_corr;

        num_pval_conf_ftest = sum(conf_ftest);
        num_pval_conf_corr  = sum(conf_corr);
        num_conf_both       = sum(conf_both);
        num_n_conf          = sum(~conf_both);
        num_adj_ftest       = sum(T1.pvalue_adj      < 0.05);
        num_adj_corr        = sum(T1.pvalue_adj_corr < 0.05);

        fprintf("  Total genes tested  : %d\n", height(T1));
        fprintf("  Confident (F+corr)  : %d\n", num_conf_both);
        fprintf("  F-test only p<0.05  : %d\n", num_pval_conf_ftest);
        fprintf("  Corr-test only p<0.05: %d\n", num_pval_conf_corr);

        % ── Sort by significance then phase then amplitude ─────────────────
        [T1, sort_idx] = sortrows(T1, ...
            ["pvalue_adj_corr","pvalue_adj","Acrophase_24","Abs_Amp"], ...
            {'ascend','ascend','ascend','descend'});
        T2 = T2(sort_idx, :);

        % ── T3: ZT00-normalised means (fallback to 2nd col if ZT00 == 0) ──
        T3       = T2;
        R0_ref   = T3{:, 2};
        zero_m   = (R0_ref == 0);
        R0_ref(zero_m) = T3{zero_m, 3};
        T3{:, 2:end} = T3{:, 2:end} ./ R0_ref;

        % ── Write ALL-genes files ──────────────────────────────────────────
        % Use the already-sanitized outdir_name as file prefix so that any
        % special characters in the cell type name (+, /, \, spaces, etc.)
        % never reach writetable and get misinterpreted as path separators.
        fbase = [outdir_name, char(per_label)];
        writetable(T1, fullfile(outdir, [fbase 'circadian_analysis_all.csv']));
        writetable(T2, fullfile(outdir, [fbase 'circadian_ZTs_mean.csv']));
        writetable(T3, fullfile(outdir, [fbase 'circadian_ZTs_mean_normalized.csv']));

        % ── Write CONFIDENT-genes files ────────────────────────────────────
        if rm_low_conf
            conf_after_sort = conf_both(sort_idx);
            T1_conf = T1(conf_after_sort, :);
            T2_conf = T2(conf_after_sort, :);
            T3_conf = T3(conf_after_sort, :);
            writetable(T1_conf, fullfile(outdir, [fbase 'circadian_analysis_confident.csv']));
            writetable(T2_conf, fullfile(outdir, [fbase 'circadian_ZTs_mean_confident.csv']));
            writetable(T3_conf, fullfile(outdir, [fbase 'circadian_ZTs_mean_normalized_confident.csv']));
        end

        % Summary row
        info_p_type(icell_type, :) = [n_cells_ct, length(sce.g), ...
                                       height(T1), num_conf_both, num_n_conf, ...
                                       num_pval_conf_ftest, num_pval_conf_corr, ...
                                       num_adj_ftest, num_adj_corr];

        if plot_heat
            generateHeatmap_circ_simple(cell_type, true, "", false, period12, outdir);
        end

    end  % cell-type loop

    % ── Cross-cell-type summary ────────────────────────────────────────────
    T0 = table(cell_type_list, ...
               info_p_type(:,1), info_p_type(:,2), info_p_type(:,3), ...
               info_p_type(:,4), info_p_type(:,5), info_p_type(:,6), ...
               info_p_type(:,7), info_p_type(:,8), info_p_type(:,9));
    T0.Properties.VariableNames = ["CellType","NumCells","NumGenes", ...
                                    "NumTested","NumConfident_Both","NumNonConfident", ...
                                    "NumPvalConf_Ftest","NumPvalConf_Corr", ...
                                    "NumAdjConf_Ftest","NumAdjConf_Corr"];

    if isempty(custom_celltype)
        if ncell_types == 1
            custom_celltype = cell_type_list{1};
        else
            custom_celltype = "all_cell_types";
        end
    end
    % Sanitise before building the file path — same rules as outdir_name so
    % special characters (+, /, \, spaces) never reach the filesystem.
    ct_safe = regexprep(strtrim(char(custom_celltype)), '[^a-zA-Z0-9_]', '_');
    ct_safe = regexprep(ct_safe, '_+', '_');
    ct_safe = regexprep(ct_safe, '^_|_$', '');
    writetable(T0, strcat(ct_safe, per_label, "summary_results.csv"));
    fprintf('\n=== TimeSCape complete.  Total time: %.1f s ===\n', toc);

end  % sce_circ_phase_estimation_stattest

% ── Local helper: Benjamini-Hochberg p-value adjustment ───────────────────
function p_adj = bh_adjust_pvalues(p)
    [p_sorted, sort_idx] = sort(p(:));
    m                    = numel(p_sorted);
    ranks                = (1 : m)';
    p_adj_sorted         = min(1, p_sorted .* m ./ ranks);
    % Enforce monotonicity (right-to-left minimum)
    for k = m-1 : -1 : 1
        p_adj_sorted(k) = min(p_adj_sorted(k), p_adj_sorted(k+1));
    end
    p_adj            = zeros(size(p));
    p_adj(sort_idx)  = p_adj_sorted;
end