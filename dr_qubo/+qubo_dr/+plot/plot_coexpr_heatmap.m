function plot_coexpr_heatmap(mat, gene_names, title_str, cmap)
% PLOT_COEXPR_HEATMAP Visualize co-expression matrix as heatmap.
%
%   PLOT_COEXPR_HEATMAP(mat, gene_names, title_str, cmap)
%
% Creates a publication-quality heatmap visualization of gene co-expression
% or other square matrices with gene names on axes.
%
% INPUT:
%   mat         - Square matrix to visualize (ng × ng)
%   gene_names  - Gene names (string array, length ng)
%   title_str   - Title for the plot (string, default '')
%   cmap        - Colormap name (string, default 'parula')
%
% AUTHOR: Selim Romero, Texas A&M University

    if nargin < 3
        title_str = '';
    end
    if nargin < 4
        cmap = 'parula';
    end

    ng = size(mat, 1);
    gene_names = string(gene_names(:))';

    % Create heatmap
    imagesc(mat);
    colorbar;
    colormap(cmap);

    % Set axis ticks and labels
    xticks(1:ng);
    yticks(1:ng);
    xticklabels(gene_names);
    yticklabels(gene_names);

    % Rotate x labels for readability
    xtickangle(90);

    % Format axes
    xlabel('Gene', 'FontWeight', 'bold');
    ylabel('Gene', 'FontWeight', 'bold');
    set(gca, 'FontSize', 8);
    set(gca, 'TickLabelInterpreter', 'none');

    % Add title
    if ~isempty(title_str)
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Ensure tight layout
    axis square;

end
