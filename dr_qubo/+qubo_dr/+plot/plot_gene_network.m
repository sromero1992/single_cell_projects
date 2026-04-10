function [h, G] = plot_gene_network(sub_Q_net, gene_names, edge_pct, test_label, reference_label)
% PLOT_GENE_NETWORK Visualize gene co-expression network with force layout.
%
%   [h, G] = PLOT_GENE_NETWORK(sub_Q_net, gene_names, edge_pct)
%   [h, G] = PLOT_GENE_NETWORK(sub_Q_net, gene_names, edge_pct, test_label, reference_label)
%
% Creates a network graph visualization of genes with edges colored by
% direction of co-expression change:
%
%   RED  — test / disease co-expression gain  (sub_Q_net(i,j) > 0)
%   BLUE — reference / control co-expression gain  (sub_Q_net(i,j) < 0)
%
% NOTE: sub_Q_net should be passed as −dS (i.e. test − reference), so
%       positive entries represent test/disease gain — the biologically
%       natural sign convention.  run_pipeline() stores this as sub_Q_net
%       automatically (sub_Q_net = −sub_Q0 = −dS).
%
% INPUT:
%   sub_Q_net      - Co-expression matrix (ng × ng), symmetric.
%                    Pass as −dS so positive = test/disease gain.
%   gene_names     - Gene names (string array, length ng)
%   edge_pct       - Percentile threshold for edge visibility (0–100, default 95)
%   test_label     - Display label for the test condition (default 'test')
%   reference_label - Display label for the reference condition (default 'reference')
%
% OUTPUT:
%   h              - Plot handle
%   G              - graph object
%
% AUTHOR: Selim Romero, Texas A&M University

    if nargin < 3; edge_pct        = 95;          end
    if nargin < 4; test_label      = 'test';       end
    if nargin < 5; reference_label = 'reference';  end

    ng = size(sub_Q_net, 1);
    gene_names = string(gene_names(:))';

    % Filter edges by percentile threshold
    edge_threshold = prctile(abs(sub_Q_net(triu(true(ng), 1))), edge_pct);

    % Create adjacency matrix with filtered edges (upper triangle only)
    adj_filtered = abs(sub_Q_net) > edge_threshold;
    adj_filtered = triu(adj_filtered, 1);

    % Build graph from upper triangle
    [src, dst] = find(adj_filtered);
    weights = sub_Q_net(sub2ind([ng, ng], src, dst));

    G = graph(src, dst, weights, ng);
    G.Nodes.Name = gene_names';

    % Create base plot with force-directed layout
    h = plot(G, 'Layout', 'force', 'Iterations', 500, ...
        'NodeColor', [0.7 0.7 0.7], 'NodeFontSize', 15, ...
        'NodeFontColor', 'black', 'MarkerSize', 6, ...
        'EdgeColor', 'k', 'LineWidth', 1);

    % Color edges by direction of co-expression change
    %   sub_Q_net > 0  →  test/disease gain   →  RED
    %   sub_Q_net < 0  →  reference gain       →  BLUE
    hold on;
    p = h.Parent;

    for i = 1:height(G.Edges)
        w = G.Edges.Weight(i);
        if w > 0
            h.EdgeCData(i) = 2;  % Red  — test/disease co-expression gain
        else
            h.EdgeCData(i) = 1;  % Blue — reference co-expression gain
        end
    end

    % Colormap: index 1 = blue (reference gain), index 2 = red (test gain)
    colormap(p, [0 0 1; 1 0 0]);
    hold off;

    % Format axes
    axis off;
    set(gcf, 'Color', 'white');

    % Title
    title(sprintf('Gene Network (%d genes, %d edges)  red=%s gain, blue=%s gain', ...
        numnodes(G), numedges(G), test_label, reference_label), ...
        'FontSize', 12, 'FontWeight', 'bold');

end
