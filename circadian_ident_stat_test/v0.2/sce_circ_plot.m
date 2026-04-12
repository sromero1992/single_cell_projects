function sce_circ_plot(sce, tmeta, cust_cells, plot_type, period12, norm_str, outdir) %#ok<INUSL>
% SCE_CIRC_PLOT
%   Generates per-gene expression plots for one cell type.
%   Reads pre-computed results from CSV files — does NOT re-run the analysis
%   or write any spreadsheets.  Run the analysis first (Step ③ in GUI).
%
% INPUTS:
%   sce        - SingleCellExperiment object (used only for API compatibility).
%   tmeta      - Tmeta table (old_labels | new_labels | ZT_times).
%   cust_cells - Cell type string.
%   plot_type  - 1 = confident genes (F-test AND corr p < 0.05)
%                2 = non-confident genes
%                3 = classical circadian gene hits
%   period12   - Logical; 12-hr (true) or 24-hr (false).
%   norm_str   - (unused, kept for API compatibility).
%
% USAGE:
%   old_labels = ["1La","2La","3La","4La","5Da","6Da","7Da","8Da"]';
%   new_labels = ["ZT00","ZT03","ZT06","ZT09","ZT12","ZT15","ZT18","ZT21"]';
%   times      = (0:3:21)';
%   tmeta      = table(old_labels, new_labels, times, 'VariableNames',{'old_labels','new_labels','ZT_times'});
%   sce_circ_plot(sce, tmeta, 'T cells', 1);

    if nargin < 4 || isempty(plot_type);  plot_type = 1;          end
    if nargin < 5 || isempty(period12);   period12  = false;      end
    if nargin < 7 || isempty(outdir)
        outdir_name = regexprep(strtrim(char(cust_cells)), '[^\w]', '_');
        outdir      = fullfile(pwd, outdir_name);
    end

    if period12; per_label = "_period_12_"; else; per_label = "_period_24_"; end

    % ── Read pre-computed results from CSV ─────────────────────────────────
    fname_stats = fullfile(outdir, sprintf('%s%scircadian_analysis_all.csv', cust_cells, per_label));
    fname_zts   = fullfile(outdir, sprintf('%s%scircadian_ZTs_mean.csv',     cust_cells, per_label));

    if ~exist(fname_stats,'file') || ~exist(fname_zts,'file')
        error(['Analysis files not found for cell type "%s".\n' ...
               'Please run the analysis first (Step ③ in the GUI).\n' ...
               'Expected: %s'], cust_cells, fname_stats);
    end

    T1 = readtable(fname_stats, 'ReadVariableNames', true);
    T2 = readtable(fname_zts,   'ReadVariableNames', true);

    fprintf("  Genes loaded from CSV: %d\n", height(T1));

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

    % ── Classical circadian gene patterns ──────────────────────────────────
    classic_circ = ["arnt","bhlh","clock","cry","dbp","tef","hlf","raf","erk", ...
                    "mek","ras","mtor","map","ral","akt","hif","kras","myc", ...
                    "nfkb","per","wnt","nr1d","rev","pik","ror"];

    ngenes     = height(T1);
    is_classic = false(ngenes, 1);
    for ig = 1:ngenes
        if any(startsWith(lower(T1.Genes(ig)), classic_circ))
            is_classic(ig) = true;
        end
    end

    if any(is_classic)
        disp("  Classical circadian genes found:");
        disp(T1.Genes(is_classic));
    else
        disp("  No classical circadian genes identified.");
    end

    % ── Gene index sets ────────────────────────────────────────────────────
    is_confident = (T1.pvalue <= 0.05) & (T1.pvalue_corr <= 0.05);
    gidx_conf    = find(is_confident)';
    gidx_nonconf = find(~is_confident)';
    gidx_classic = find(is_classic)';

    % ── Output paths (subdirectories inside the cell-type folder) ────────
    base    = char(strcat(cust_cells, per_label));
    p_conf  = fullfile(outdir, [base 'plots_confident']);
    p_nconf = fullfile(outdir, [base 'plots_non_confident']);
    p_circ  = fullfile(outdir, [base 'plots_classic_circadian']);

    % ── Helper: save plot for one gene ─────────────────────────────────────
    function save_gene_plot(i, out_path)
        % Force white theme regardless of MATLAB colour scheme
        f  = figure('Visible','off','Color',[1 1 1],'InvertHardcopy','off');
        ax = axes('Parent',f);

        % ── White axes, black labels everywhere ──────────────────────────
        set(ax, ...
            'Color',         [1 1 1], ...
            'XColor',        [0 0 0], ...
            'YColor',        [0 0 0], ...
            'GridColor',     [0.85 0.85 0.85], ...
            'GridAlpha',     0.6, ...
            'MinorGridColor',[0.90 0.90 0.90], ...
            'FontSize',      11, ...
            'FontName',      'Helvetica', ...
            'Box',           'on', ...
            'LineWidth',     0.8);

        amp_i  = T1.Amp(i);
        acro_i = T1.Acrophase(i);
        per_i  = T1.Period(i);
        mes_i  = T1.Mesor(i);
        fval_i = amp_i * cos(2*pi*(tval - acro_i) / per_i) + mes_i;
        Rzts_i = table2array(T2(i, 2:end));

        plot(ax, tval, fval_i, '-',  'Color',[0.13 0.47 0.71], 'LineWidth', 2);
        hold(ax, 'on');
        plot(ax, actual_times, Rzts_i, 'o-', ...
             'Color',[0.85 0.33 0.10], 'LineWidth', 1.8, 'MarkerSize', 7, ...
             'MarkerFaceColor',[0.85 0.33 0.10]);
        xl = xline(ax, T1.Acrophase_24(i), '--', 'Color',[0.6 0 0], 'LineWidth', 1.2);
        xl.Label = sprintf('ZT%.1f', T1.Acrophase_24(i));
        xl.LabelVerticalAlignment = 'bottom';
        xl.FontSize = 9;
        grid(ax, 'on');
        ylim(ax, [0 inf]);
        xlim(ax, [t0 tf]);
        xticks(ax, actual_times);
        xticklabels(ax, string(actual_times));

        title(ax, sprintf('%s  |  F-test p=%.3g  |  Corr p=%.3g', ...
              T1.Genes{i}, T1.pvalue(i), T1.pvalue_corr(i)), ...
              'Color','k', 'FontSize', 12, 'FontWeight','bold');
        xlabel(ax, 'Zeitgeber Time (hrs)',       'Color','k', 'FontSize', 11);
        ylabel(ax, 'Expression (log-normalised)','Color','k', 'FontSize', 11);

        leg = legend(ax, {'Cosine fit','Mean / ZT','Acrophase'}, ...
                     'Location','northwest','FontSize',10);
        set(leg, 'TextColor','k','EdgeColor',[0.7 0.7 0.7],'Color',[1 1 1]);

        hold(ax,'off');
        fname_png = fullfile(out_path, strcat("plot_", T1.Genes(i), ".png"));
        print(f, fname_png, '-dpng', '-r150');   % explicit white-bg PNG export
        close(f);
    end

    % ── Plot type 1: Confident genes ───────────────────────────────────────
    if plot_type == 1
        mkdir(p_conf);
        fprintf("  Saving %d confident-gene plots → %s\n", length(gidx_conf), p_conf);
        for i = gidx_conf; save_gene_plot(i, p_conf); end
    end

    % ── Plot type 2: Non-confident genes ──────────────────────────────────
    if plot_type == 2
        mkdir(p_nconf);
        fprintf("  Saving %d non-confident-gene plots → %s\n", length(gidx_nonconf), p_nconf);
        for i = gidx_nonconf; save_gene_plot(i, p_nconf); end
    end

    % ── Plot type 3: Classical circadian genes ─────────────────────────────
    if plot_type == 3
        mkdir(p_circ);
        fprintf("  Saving %d classic-circadian-gene plots → %s\n", length(gidx_classic), p_circ);
        for i = gidx_classic; save_gene_plot(i, p_circ); end
    end

end
