function scReservoir_plot_landscape(pc1, pc2, energy, prob, pseudotime, attractorCells)
% SCRESERVOIR_PLOT_LANDSCAPE  Plot Waddington energy landscape and fate probabilities.
%
%   scReservoir_plot_landscape(pc1, pc2, energy, prob, pseudotime, attractorCells)
%
%   Called automatically by scReservoir_landscape when 'plot' is true.
%   Can also be called manually for custom visualization.

nAttractors = size(prob, 2);

figure('Name', 'scReservoir: Waddington Landscape', 'NumberTitle', 'off', ...
       'Position', [100 100 1400 500]);

% --- Panel 1: Energy landscape (3D) -------------------------------------
subplot(1, 3, 1);
scatter3(pc1, pc2, energy, 12, pseudotime, 'filled', 'MarkerFaceAlpha', 0.6);
colormap(gca, 'viridis');
colorbar;
hold on;
% Mark attractors
scatter3(pc1(attractorCells), pc2(attractorCells), energy(attractorCells), ...
         60, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
hold off;
xlabel('PC1'); ylabel('PC2'); zlabel('Energy (−log P)');
title('Waddington Energy Landscape');
legend({'Cells (pseudotime)', 'Attractors'}, 'Location', 'best');
grid on; view([-30 30]);

% --- Panel 2: Cell density / potential (2D heatmap) --------------------
subplot(1, 3, 2);
% Create grid interpolation of energy
xi = linspace(min(pc1), max(pc1), 80);
yi = linspace(min(pc2), max(pc2), 80);
[XI, YI] = meshgrid(xi, yi);
try
    EI = griddata(pc1, pc2, energy, XI, YI, 'natural');
    contourf(XI, YI, EI, 20, 'LineColor', 'none');
    colormap(gca, flipud(hot));
    colorbar;
    hold on;
    scatter(pc1(attractorCells), pc2(attractorCells), 50, 'w', ...
            'filled', 'MarkerEdgeColor', 'k');
    hold off;
catch
    scatter(pc1, pc2, 10, energy, 'filled');
    colormap(gca, flipud(hot)); colorbar;
end
xlabel('PC1'); ylabel('PC2');
title('Energy Potential (top view)');

% --- Panel 3: Fate probabilities (one per attractor) --------------------
subplot(1, 3, 3);
% Show probability of most dominant fate per cell
[maxProb, dominantFate] = max(prob, [], 2);
scatter(pc1, pc2, 15, dominantFate, 'filled', 'MarkerFaceAlpha', 0.7);
colormap(gca, lines(nAttractors));
caxis([1 nAttractors]);
cb = colorbar;
cb.Ticks = 1:nAttractors;
cb.TickLabels = arrayfun(@(i) sprintf('Fate %d', i), 1:nAttractors, 'UniformOutput', false);
xlabel('PC1'); ylabel('PC2');
title('Dominant Cell Fate');

sgtitle('scReservoir: Dynamical Landscape Analysis', 'FontSize', 14, 'FontWeight', 'bold');
end
