function generateHeatmap_circ_simple(celltype, strict, customName, circ, period12, outdir, axHandle)
% GENERATEHEATMAP_CIRC_SIMPLE
%   Reads circadian analysis CSV files, filters and sorts genes, and renders
%   a blue-white-red z-score heatmap (genes × ZT time points).
%   Can display in a supplied GUI axes AND/OR save to the cell-type directory.
%
% INPUTS:
%   celltype   - String; cell-type label matching the analysis file names.
%   strict     - Logical; true = keep only genes with pvalue AND pvalue_corr < 0.05.
%   customName - String appended to saved file names ('' for none).
%   circ       - Logical; true = restrict to classical circadian gene prefixes.
%   period12   - (default false) Logical; true = read 12-hr files.
%   outdir     - (default pwd) Directory where CSV files were written and
%                where the heatmap PNG/FIG will be saved.
%   axHandle   - (default []) GUI axes handle.  When supplied the heatmap is
%                also rendered live into the GUI panel.

    if nargin < 5 || isempty(period12); period12  = false; end
    if nargin < 6 || isempty(outdir);   outdir    = pwd;   end
    if nargin < 7;                      axHandle  = [];    end

    % ── Locate input files ─────────────────────────────────────────────────
    if period12; suffix = '_period_12_'; else; suffix = '_period_24_'; end
    fname  = fullfile(outdir, sprintf('%s%scircadian_analysis_all.csv', celltype, suffix));
    fname2 = fullfile(outdir, sprintf('%s%scircadian_ZTs_mean.csv',     celltype, suffix));

    if ~exist(fname, 'file') || ~exist(fname2, 'file')
        warning('Analysis file not found:\n  %s\nRun the analysis first.', fname);
        return;
    end

    D    = readtable(fname,  'ReadVariableNames', true);
    Dzts = readtable(fname2, 'ReadVariableNames', true);

    % ── Sort by acrophase → amplitude ─────────────────────────────────────
    [D, idx] = sortrows(D, ["Acrophase_24","Abs_Amp"], ["ascend","descend"]);
    Dzts     = Dzts(idx, :);

    % ── Filter to confident genes ──────────────────────────────────────────
    if strict
        keep = (D.pvalue < 0.05) & (D.pvalue_corr < 0.05);
        D    = D(keep, :);
        Dzts = Dzts(keep, :);
    end

    % ── Correct negative-amplitude rows ───────────────────────────────────
    neg_amp                  = D.Amp < 0;
    D.Amp(neg_amp)           = abs(D.Amp(neg_amp));
    D.Acrophase(neg_amp)     = D.Acrophase(neg_amp) + 12;
    D.Acrophase_24(neg_amp)  = mod(D.Acrophase_24(neg_amp) + 12, 24);

    [D, idx2] = sortrows(D, ["Acrophase_24","Abs_Amp"], ["ascend","descend"]);
    Dzts      = Dzts(idx2, :);

    % ── Restrict to classical circadian genes ─────────────────────────────
    if circ
        classic_circ = ["arnt","bhlh","clock","cry","dbp","tef","hlf","raf","erk", ...
                        "mek","ras","mtor","map","ral","akt","hif","kras","myc", ...
                        "nfkb","per","wnt","nr1d","rev","pik","ror"];
        glcirc = contains(lower(D.Genes), classic_circ, 'IgnoreCase', true);
        D    = D(glcirc, :);
        Dzts = Dzts(glcirc, :);
    end

    if height(D) <= 1
        fprintf('Too few genes after filtering (%d) – skipping heatmap.\n', height(D));
        return;
    end

    % ── Build expression matrix ────────────────────────────────────────────
    zt_cols      = 2 : width(Dzts);
    time_batches = string(Dzts.Properties.VariableNames(zt_cols));
    gzts         = Dzts{:, zt_cols};
    data_scaled  = normalize(gzts, 2, 'zscore');

    % ── Colormap (blue→white→red) ──────────────────────────────────────────
    n           = 256;
    blueToWhite = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
    whiteToRed  = [ones(n/2,1),        linspace(1,0,n/2)', linspace(1,0,n/2)'];
    cmap        = [blueToWhite; whiteToRed];

    % ── Shared render function ─────────────────────────────────────────────
    function paint_into(ax, parent_fig)
        cla(ax);
        h = imagesc(ax, data_scaled);
        set(h, 'AlphaData', ~isnan(data_scaled));
        set(ax, 'Color','white', 'XColor',[0.12 0.12 0.12], 'YColor',[0.12 0.12 0.12]);
        ax.XTick       = 1 : length(time_batches);
        ax.XTickLabels = time_batches;
        ax.XTickLabelRotation = 45;
        if height(D) <= 60
            ax.YTick       = 1 : height(D);
            ax.YTickLabels = string(Dzts{:,1});
            ax.FontSize    = max(6, min(10, floor(300 / height(D))));
        else
            ax.YTick       = [];
            ax.YTickLabels = {};
            ax.FontSize    = 9;
        end
        title(ax, sprintf('%s  –  %d genes  (z-score)', celltype, height(D)), ...
              'Color','k', 'FontWeight','bold', 'FontSize', 11);
        xlabel(ax, 'Zeitgeber Time', 'Color','k', 'FontSize', 10);
        ylabel(ax, 'Genes  (sorted by acrophase)', 'Color','k', 'FontSize', 10);
        dr = [min(data_scaled(:)), max(data_scaled(:))];
        if diff(dr) > 0; clim(ax, dr); end
        colormap(parent_fig, cmap);
        colorbar(ax);
        axis(ax, 'tight');
    end

    % ── Render into GUI axes (live preview) ───────────────────────────────
    if ~isempty(axHandle) && ishandle(axHandle)
        paint_into(axHandle, ancestor(axHandle, 'figure'));
        drawnow;
    end

    % ── Build output file name prefix ─────────────────────────────────────
    fp = char(celltype);
    if ~isempty(customName) && ~strcmp(customName,''); fp = [fp '_' char(customName)]; end
    fp = [fp '_'];
    if strict; fp = [fp 'strict']; else; fp = [fp 'all']; end
    if circ;   fp = [fp '_core_circ']; end

    % ── Save CSV of filtered/sorted gene table ────────────────────────────
    writetable(D, fullfile(outdir, [fp '_heatmap_genes.csv']));

    % ── Save PNG + FIG ────────────────────────────────────────────────────
    figSave = figure('Visible','off','Color','white','Position',[0 0 900 600]);
    ax_save = axes(figSave);
    paint_into(ax_save, figSave);

    out_png = fullfile(outdir, [fp '_heatmap.png']);
    out_fig = fullfile(outdir, [fp '_heatmap.fig']);
    saveas(figSave, out_png);
    saveas(figSave, out_fig);
    close(figSave);
    fprintf('Heatmap saved → %s\n', out_png);

end  % generateHeatmap_circ_simple
