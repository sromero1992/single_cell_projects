function sce_circ_plot_gene(sce, tmeta, cust_cells, period12, cust_gene, ...
                            axHandle, print_scdata, norm_str, use_violin_plot, outdir)
% SCE_CIRC_PLOT_GENE
%   Plots the cosine fit and mean-per-timepoint expression for one gene
%   in one cell type.  Reads pre-computed stats from CSV — does NOT
%   re-run the analysis or write any spreadsheets.
%   Run the analysis first (Step ③ in the GUI).
%
% INPUTS:
%   sce            - SingleCellExperiment object (needed only when
%                    print_scdata = true for single-cell overlay).
%   tmeta          - Tmeta table (old_labels | new_labels | ZT_times).
%   cust_cells     - Cell type string to restrict to.
%   period12       - Logical; 12-hr (true) or 24-hr (false).
%   cust_gene      - Gene name string.
%   axHandle       - Axes handle to plot into.
%   print_scdata   - Logical; overlay single-cell scatter/violin.
%   norm_str       - 'lib_size' | 'magic_impute' | 'none'.
%                    'none' = use sce.X as-is (no transformation).
%   use_violin_plot- Logical; violin (true) or scatter (false).
%   outdir         - Full path to cell-type output directory.

    if nargin < 4  || isempty(period12);        period12        = false;     end
    if nargin < 5  || isempty(cust_gene);       error('Gene name required.'); end
    if nargin < 6  || isempty(axHandle);        error('Axes handle required.'); end
    if nargin < 7  || isempty(print_scdata);    print_scdata    = false;     end
    if nargin < 8  || isempty(norm_str);        norm_str        = 'lib_size'; end
    if nargin < 9  || isempty(use_violin_plot); use_violin_plot = false;     end
    if nargin < 10 || isempty(outdir)
        outdir_name = regexprep(strtrim(char(cust_cells)), '[^\w]', '_');
        outdir      = fullfile(pwd, outdir_name);
    end

    if period12; period = 12; per_label = '_period_12_';
    else;        period = 24; per_label = '_period_24_'; end

    % ── Read pre-computed results from CSV ─────────────────────────────────
    fname_stats = fullfile(outdir, sprintf('%s%scircadian_analysis_all.csv', cust_cells, per_label));
    fname_zts   = fullfile(outdir, sprintf('%s%scircadian_ZTs_mean.csv',     cust_cells, per_label));

    if ~exist(fname_stats,'file') || ~exist(fname_zts,'file')
        error(['Analysis files not found for cell type "%s" (period=%d hr).\n' ...
               'Please run the analysis first (Step ③ in the GUI).\n' ...
               'Expected: %s'], cust_cells, period, fname_stats);
    end

    T1 = readtable(fname_stats, 'ReadVariableNames', true);
    T2 = readtable(fname_zts,   'ReadVariableNames', true);

    gene_idx = find(strcmp(T1.Genes, cust_gene));
    if isempty(gene_idx)
        error('Gene "%s" not found in results for cell type "%s".\nCheck spelling or run the analysis first.', ...
              cust_gene, cust_cells);
    end

    % ── Identify actual time points from T2 column names ──────────────────
    zt_labels    = T2.Properties.VariableNames(2:end);
    actual_times = zeros(1, length(zt_labels));
    for it = 1:length(zt_labels)
        idx_t = find(tmeta.new_labels == string(zt_labels{it}), 1);
        if ~isempty(idx_t); actual_times(it) = tmeta.ZT_times(idx_t); end
    end
    t0   = min(actual_times);
    tf   = max(actual_times);
    tval = t0 : 0.1 : tf;

    % ── Fitted cosine curve ────────────────────────────────────────────────
    amp_g   = T1.Amp(gene_idx);
    acro_g  = T1.Acrophase(gene_idx);
    mesor_g = T1.Mesor(gene_idx);
    fval    = amp_g * cos(2*pi*(tval - acro_g) / period) + mesor_g;

    Rzts = table2array(T2(gene_idx, 2:end));

    % ── Always count total cells for N in title (no normalization needed) ──
    n_cells_type = sum(sce.c_cell_type_tx == cust_cells);

    % ── Collect per-cell expression (only when overlay requested) ─────────
    nc_cum   = 0;
    all_expr = [];
    all_tp   = [];

    if print_scdata
        if strcmp(norm_str, 'lib_size')
            X_norm = pkg.norm_libsize(sce.X, 1e4);
            X_norm = log1p(X_norm);
            sce.X  = X_norm; clear X_norm;
        elseif strcmp(norm_str, 'magic_impute')
            X_norm = sc_impute(sce.X, 'MAGIC');
            sce.X  = X_norm; clear X_norm;
        end
        % 'none': sce.X already holds data as-is — no transformation applied.

        ic0     = find(sce.c_cell_type_tx == cust_cells);
        sce_sub = sce.selectcells(ic0);
        ig      = find(sce_sub.g == cust_gene);
        clear ic0;

        batch_time = unique(sce_sub.c_batch_id);
        nzts_sc    = length(batch_time);
        at_sc      = nan(nzts_sc, 1);
        for it = 1:nzts_sc
            idx_t = find(tmeta.new_labels == batch_time(it), 1);
            if ~isempty(idx_t); at_sc(it) = tmeta.ZT_times(idx_t); end
        end
        valid_sc   = ~isnan(at_sc);
        batch_time = batch_time(valid_sc);
        at_sc      = at_sc(valid_sc);
        nzts_sc    = sum(valid_sc);

        ncell_total = size(sce_sub.X, 2);
        all_expr    = zeros(ncell_total, 1);
        all_tp      = zeros(ncell_total, 1);
        for it = 1:nzts_sc
            ics    = find(sce_sub.c_batch_id == batch_time(it));
            nc_loc = length(ics);
            if nc_loc > 0 && ~isempty(ig)
                all_expr(nc_cum+1 : nc_cum+nc_loc) = full(sce_sub.X(ig, ics));
                all_tp  (nc_cum+1 : nc_cum+nc_loc) = at_sc(it);
            end
            nc_cum = nc_cum + nc_loc;
        end
        all_expr = all_expr(1:nc_cum);
        all_tp   = all_tp  (1:nc_cum);
        clear sce_sub;
    end

    % ── Helper: apply white axes theme ────────────────────────────────────
    function apply_white_theme()
        set(ancestor(axHandle,'figure'), 'Color',[1 1 1]);
        set(axHandle, ...
            'Color',     [1 1 1], ...
            'XColor',    [0 0 0], ...
            'YColor',    [0 0 0], ...
            'GridColor', [0.82 0.82 0.82], ...
            'GridAlpha', 0.7, ...
            'FontSize',  11, ...
            'FontName',  'Helvetica', ...
            'Box',       'on');
    end

    % ── Prepare axes ───────────────────────────────────────────────────────
    apply_white_theme();
    axes(axHandle);
    cla(axHandle,'reset');
    apply_white_theme();

    % Explicit legend handle accumulators
    leg_handles = gobjects(0);
    leg_labels  = {};

    if print_scdata && nc_cum > 0 && use_violin_plot
        % ── VIOLIN MODE ────────────────────────────────────────────────────
        % Clip y-axis so extreme outliers don't flatten the cosine/mean view
        y_clip = max(prctile(all_expr, 99) * 1.15, max(Rzts) * 3);

        tbl_violin = table(all_tp, all_expr, ...
                           'VariableNames', {'TimePoints','Expression'});
        violinplot(tbl_violin, "TimePoints", "Expression");

        % violinplot overrides axes colors — restore white immediately, then
        % hide ALL violin children from the legend so our own handles dominate
        apply_white_theme();
        set(axHandle.Children, 'HandleVisibility', 'off');
        hold(axHandle, 'on');

        % Cosine fit — on top with thick line
        h_cos = plot(axHandle, tval, fval, '-', ...
                     'Color',[0.13 0.47 0.71], 'LineWidth', 3);

        % Mean per ZT — prominent red markers
        h_mean = plot(axHandle, actual_times, Rzts, 'o-', ...
                      'Color',[0.85 0.10 0.10], 'LineWidth', 2.5, ...
                      'MarkerSize', 10, ...
                      'MarkerFaceColor',[0.85 0.10 0.10], ...
                      'MarkerEdgeColor','w');

        % Acrophase vertical line — hidden from legend (label printed on axis)
        xl = xline(axHandle, T1.Acrophase_24(gene_idx), '--', ...
                   'Color',[0.5 0 0.5], 'LineWidth', 1.8, ...
                   'HandleVisibility','off');
        xl.Label = sprintf('ZT%.1f', T1.Acrophase_24(gene_idx));
        xl.LabelVerticalAlignment = 'bottom';
        xl.FontSize = 9;

        % Dummy handles for legend swatches (plotted at NaN so invisible)
        h_violin = patch(axHandle, NaN, NaN, [0.60 0.81 0.93], ...
                         'FaceAlpha', 0.55, 'EdgeColor',[0.13 0.47 0.71], ...
                         'LineWidth', 1.2);
        h_acro_v = plot(axHandle, NaN, NaN, '--', ...
                        'Color',[0.5 0 0.5], 'LineWidth', 1.8);

        grid(axHandle, 'on');
        ylim(axHandle, [0, y_clip]);
        xlim(axHandle, [t0 tf]);
        xticks(axHandle, actual_times);
        xticklabels(axHandle, string(actual_times));

        % Legend: Violin swatch first, then cosine, mean, acrophase
        leg_handles = [h_violin, h_cos, h_mean, h_acro_v];
        leg_labels  = {'Violin', 'Cosine fit', 'Mean per ZT', 'Acrophase'};

    else
        % ── NORMAL / DOTS MODE ─────────────────────────────────────────────
        hold(axHandle, 'on');

        % Cosine fit
        h_cos = plot(axHandle, tval, fval, '-', ...
                     'Color',[0.13 0.47 0.71], 'LineWidth', 2);

        % Mean per ZT
        h_mean = plot(axHandle, actual_times, Rzts, 'o-', ...
                      'Color',[0.85 0.33 0.10], 'LineWidth', 1.8, ...
                      'MarkerSize', 7, ...
                      'MarkerFaceColor',[0.85 0.33 0.10]);

        leg_handles = [h_cos, h_mean];
        leg_labels  = {'Cosine fit', 'Mean per ZT'};

        if print_scdata && nc_cum > 0
            % Single-cell scatter overlay
            h_sc = scatter(axHandle, all_tp, all_expr, 18, [0.3 0.5 0.9], ...
                           'filled', 'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22);
            leg_handles(end+1) = h_sc;
            leg_labels{end+1}  = 'Single cells';
        end

        % Acrophase vertical line — hidden from legend (label printed on axis)
        xl = xline(axHandle, T1.Acrophase_24(gene_idx), '--', ...
                   'Color',[0.6 0 0], 'LineWidth', 1.2, ...
                   'HandleVisibility','off');
        xl.Label = sprintf('ZT%.1f', T1.Acrophase_24(gene_idx));
        xl.LabelVerticalAlignment = 'bottom';
        xl.FontSize = 9;

        % Dummy handle for acrophase legend swatch
        h_acro_d = plot(axHandle, NaN, NaN, '--', ...
                        'Color',[0.6 0 0], 'LineWidth', 1.2);
        leg_handles(end+1) = h_acro_d;
        leg_labels{end+1}  = 'Acrophase';

        grid(axHandle, 'on');
        ylim(axHandle, [0 inf]);
        xlim(axHandle, [t0 tf]);
        xticks(axHandle, actual_times);
        xticklabels(axHandle, string(actual_times));
    end

    % ── Title, labels, legend ──────────────────────────────────────────────
    title_str = sprintf('%s  |  F-test p=%.3g  |  Corr p=%.3g  |  Phase ZT%.1f  |  N=%d', ...
        T1.Genes{gene_idx}, T1.pvalue(gene_idx), T1.pvalue_corr(gene_idx), ...
        T1.Acrophase_24(gene_idx), n_cells_type);
    title(axHandle, title_str, 'Color','k', 'FontSize', 12, 'FontWeight','bold');
    xlabel(axHandle, 'Zeitgeber Time (hrs)',        'Color','k', 'FontSize', 11);
    ylabel(axHandle, 'Expression (log-normalised)', 'Color','k', 'FontSize', 11);
    set(axHandle, 'FontSize', 11);

    % Legend with explicit object handles → labels always match the right line
    leg = legend(axHandle, leg_handles, leg_labels, ...
                 'Location','northwest', 'FontSize', 10);
    set(leg, 'TextColor','k', 'EdgeColor',[0.7 0.7 0.7], 'Color',[1 1 1]);

    hold(axHandle, 'off');

end
