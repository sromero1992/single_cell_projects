function [h, G] = plot_gene_network(sub_Q_net, gene_names, edge_pct)
% PLOT_GENE_NETWORK Visualize gene co-expression network with force layout.
%
%   [h, G] = PLOT_GENE_NETWORK(sub_Q_net, gene_names, edge_pct)
%
% Creates a network graph visualization of genes with edges colored by
% type (blue = WT co-expression gain, red = KO co-expression gain).
% Uses force-directed layout for intuitive visualization.
%
% INPUT:
%   sub_Q_net    - Co-expression matrix (ng × ng), symmetric
%   gene_names   - Gene names (string array, length ng)
%   edge_pct     - Percentile threshold for edge visibility (0-100, default 95)
%
% OUTPUT:
%   h            - Plot handle
%   G            - graph object
%
% AUTHOR: Selim Romero, Texas A&M University

    if nargin < 3
        edge_pct = 95;
    end

    ng = size(sub_Q_net, 1);
    gene_names = string(gene_names(:))';

    % Filter edges by percentile threshold
    edge_threshold = prctile(abs(sub_Q_net(triu(true(ng), 1))), edge_pct);

    % Create adjacency matrix with filtered edges
    adj_filtered = abs(sub_Q_net) > edge_threshold;
    adj_filtered = triu(adj_filtered, 1);

    % Build graph from upper triangle
    [src, dst] = find(adj_filtered);
    weights = sub_Q_net(sub2ind([ng, ng], src, dst));

    G = graph(src, dst, weights, ng);
    G.Nodes.Name = gene_names';

    % Create plot
    h = plot(G, 'Layout', 'force', 'Iterations', 500, ...
        'NodeColor', [0.7 0.7 0.7], 'NodeFontSize', 15, ...
        'NodeFontColor', 'black', 'MarkerSize', 6, ...
        'EdgeColor', 'k', 'LineWidth', 1);

    % Color edges by type (positive = WT, negative = KO)
    edge_data = table(G.Edges.EndNodes, G.Edges.Weight, ...
        'VariableNames', {'EndNodes', 'Weight'});

    hold on;
    p = h.Parent;
    % Clear current edges for recoloring
    edge_handles = h.EdgeCurveAlpha;

    % Recolor edges
    for i = 1:height(G.Edges)
        w = G.Edges.Weight(i);
        if w > 0
            % WT co-expression (blue)
            h.EdgeCData(i) = 1;  % Blue
        else
            % KO co-expression (red)
            h.EdgeCData(i) = 2;  % Red
        end
    end

    % Set edge colormap
    colormap(p, [0 0 1; 1 0 0]);  % Blue and Red
    hold off;

    % Format axes
    axis off;
    set(gcf, 'Color', 'white');

    % Add title with information
    title(sprintf('Gene Network (%d genes, %d edges)', ...
        numnodes(G), numedges(G)), ...
        'FontSize', 12, 'FontWeight', 'bold');

end
