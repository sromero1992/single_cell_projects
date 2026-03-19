function scReservoir_plot_GRN(GRN, geneNames, varargin)
% SCRESERVOIR_PLOT_GRN  Visualize inferred gene regulatory network.
%
%   scReservoir_plot_GRN(GRN, geneNames)
%   scReservoir_plot_GRN(GRN, geneNames, 'threshold', 0.2, 'topN', 30)
%
%   INPUTS
%     GRN        nGenes x nGenes normalized influence matrix from scReservoir_GRN
%     geneNames  1 x nGenes cell array of gene names
%
%   OPTIONAL NAME-VALUE PAIRS
%     'threshold'  Edge cutoff for binary display (default: 0.1)
%     'topN'       Show only top N genes by out-degree (default: 50)
%     'cmap'       Colormap for heatmap               (default: 'hot')
%     'show_net'   Show network graph in addition to heatmap (default: true)

p = inputParser;
addRequired(p,  'GRN',       @isnumeric);
addRequired(p,  'geneNames', @iscell);
addParameter(p, 'threshold', 0.1,   @isnumeric);
addParameter(p, 'topN',      50,    @isnumeric);
addParameter(p, 'cmap',      'hot', @ischar);
addParameter(p, 'show_net',  true,  @islogical);
parse(p, GRN, geneNames, varargin{:});
opt = p.Results;

nGenes = size(GRN, 1);
topN   = min(opt.topN, nGenes);

% Select top genes by total out-degree
out_degree = sum(GRN > opt.threshold, 2);
[~, topIdx] = sort(out_degree, 'descend');
topIdx = topIdx(1:topN);
GRN_sub   = GRN(topIdx, topIdx);
names_sub = geneNames(topIdx);

% --- Figure setup --------------------------------------------------------
if opt.show_net
    figure('Name','scReservoir: GRN','NumberTitle','off','Position',[100 100 1400 600]);
    nPanels = 2;
else
    figure('Name','scReservoir: GRN','NumberTitle','off','Position',[100 100 700 600]);
    nPanels = 1;
end

% --- Panel 1: Heatmap ---------------------------------------------------
subplot(1, nPanels, 1);
imagesc(GRN_sub);
colormap(gca, opt.cmap);
colorbar;
caxis([0 1]);
xlabel('Target Gene'); ylabel('Regulator Gene');
title(sprintf('GRN (top %d genes by out-degree)', topN));
set(gca, 'XTick', 1:topN, 'XTickLabel', names_sub, 'XTickLabelRotation', 45, ...
         'YTick', 1:topN, 'YTickLabel', names_sub, 'FontSize', 7);

% --- Panel 2: Network graph (if requested) ------------------------------
if opt.show_net
    subplot(1, 2, 2);
    GRN_bin = GRN_sub > opt.threshold;
    GRN_bin(logical(eye(topN))) = 0;

    % Build digraph
    [src, tgt] = find(GRN_bin);
    weights    = GRN_sub(sub2ind(size(GRN_sub), src, tgt));

    if ~isempty(src)
        G = digraph(src, tgt, weights, names_sub);
        deg = indegree(G) + outdegree(G);
        node_sizes = 3 + 15 * (deg / max(deg + 1));

        h_plot = plot(G, ...
            'Layout',    'force', ...
            'NodeLabel', names_sub, ...
            'MarkerSize', node_sizes, ...
            'ArrowSize',  8, ...
            'LineWidth',  1.5 * weights / max(weights), ...
            'NodeColor', [0.2 0.5 0.8], ...
            'EdgeColor', [0.7 0.2 0.2]);

        title(sprintf('GRN Network (threshold=%.2f)', opt.threshold));
    else
        text(0.5, 0.5, 'No edges above threshold', ...
             'HorizontalAlignment','center','Units','normalized');
        title('GRN Network (no edges)');
    end
    axis off;
end

sgtitle('scReservoir: Gene Regulatory Network', 'FontSize', 13, 'FontWeight', 'bold');
end
