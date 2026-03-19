%% example03_attractor_landscape.m
%
% scReservoir: Full Waddington Landscape Pipeline
%
% Demonstrates the complete pipeline:
%   1. Preprocess scRNA-seq data
%   2. Run sparse reservoir + randomized SVD
%   3. Estimate latent dynamics (dh/dt = A*h)
%   4. Detect cell fate attractors
%   5. Reconstruct Waddington energy landscape
%   6. Compute fate probabilities
%   7. Infer GRN from latent dynamics
%   8. Identify marker genes per attractor
%
% For atlas-scale data (>50k cells), use scReservoir_scalable directly.

clear; clc; close all;
addpath(genpath('../scReservoir'));

%% =====================================================================
%  1. Load / Simulate Multi-Fate Data
%% =====================================================================
rng(2025);
nCells  = 1000;
nGenes  = 400;
nFates  = 4;

fprintf('Simulating %d-fate differentiation dataset...\n', nFates);

% Simulate cells spreading toward nFates attractors from a common progenitor
pseudotime = rand(nCells, 1);
fate_labels = ceil(rand(nCells, 1) * nFates);

% Attractor centers in gene space
centers = randn(nFates, nGenes) * 2;

X_raw = zeros(nCells, nGenes);
for c = 1:nCells
    pt = pseudotime(c);
    f  = fate_labels(c);
    % Blend from origin toward fate attractor
    X_raw(c, :) = (1 - pt) * randn(1, nGenes) * 0.5 + pt * centers(f, :) + ...
                  randn(1, nGenes) * 0.3;
end
X_raw = max(0, X_raw);

geneNames = arrayfun(@(i) sprintf('Gene%04d', i), 1:nGenes, 'UniformOutput', false);

%% =====================================================================
%  2. Preprocess
%% =====================================================================
[X, geneNames, ~] = scReservoir_preprocess(X_raw, geneNames, ...
    'nHVG', 250, 'normalize', true, 'log_transform', true);
[nCells, nGenes] = size(X);
fprintf('Preprocessed: %d cells x %d genes\n', nCells, nGenes);

%% =====================================================================
%  3. Full Landscape Pipeline
%% =====================================================================
fprintf('\nRunning scReservoir landscape pipeline...\n');

result = scReservoir_landscape(X, pseudotime, geneNames, ...
    'N_res',         500, ...
    'rho',           0.9, ...
    'leak_rate',     0.3, ...
    'density',       0.01, ...
    'input_scale',   0.3, ...
    'rankSVD',       40, ...
    'lambda',        1e-3, ...
    'nAttractors',   nFates, ...
    'attractor_pct', 5, ...
    'k',             5, ...
    'seed',          1, ...
    'plot',          true);

%% =====================================================================
%  4. Evaluate Attractor Recovery
%% =====================================================================
fprintf('\n--- Attractor Analysis ---\n');
fprintf('Cells identified as attractors: %d (%.1f%%)\n', ...
    numel(result.attractorCells), ...
    100 * numel(result.attractorCells) / nCells);

% Check concordance with ground-truth fate labels (if available)
true_fates_at_attractors = fate_labels(result.attractorCells);
fprintf('Ground-truth fate distribution at attractors:\n');
for f = 1:nFates
    pct = 100 * sum(true_fates_at_attractors == f) / numel(result.attractorCells);
    fprintf('  Fate %d: %.1f%%\n', f, pct);
end

%% =====================================================================
%  5. Fate Probability Heatmap
%% =====================================================================
figure('Name','Fate Probabilities','Position',[100 100 900 500]);

subplot(1, 2, 1);
imagesc(result.fateProbabilities(1:100, :)');
colormap(hot); colorbar;
xlabel('Cell (first 100)'); ylabel('Attractor Fate');
title('Fate Probabilities (first 100 cells)');
yticks(1:nFates);
yticklabels(arrayfun(@(i) sprintf('Fate %d', i), 1:nFates, 'UniformOutput', false));

subplot(1, 2, 2);
[~, dom_fate] = max(result.fateProbabilities, [], 2);
pc = result.pca_scores;
gscatter(pc(:,1), pc(:,2), dom_fate, lines(nFates), '.', 12);
xlabel('PC1'); ylabel('PC2');
title('Dominant Fate Assignment in Latent Space');
legend(arrayfun(@(i) sprintf('Fate %d', i), 1:nFates, 'UniformOutput', false));

%% =====================================================================
%  6. Marker Genes Per Attractor
%% =====================================================================
fprintf('\n--- Top marker genes per attractor ---\n');
for f = 1:nFates
    % Find genes with highest expression in this attractor vs others
    mean_this  = result.attractorGenes(f, :);
    mean_other = mean(result.attractorGenes(setdiff(1:nFates, f), :), 1);
    fold_change = mean_this ./ (mean_other + 1e-8);
    [~, topMarkerIdx] = sort(fold_change, 'descend');
    topMarkerIdx = topMarkerIdx(1:5);
    fprintf('Fate %d top markers: ', f);
    fprintf('%s, ', geneNames{topMarkerIdx});
    fprintf('\n');
end

%% =====================================================================
%  7. GRN Summary
%% =====================================================================
fprintf('\n--- GRN from Latent Dynamics ---\n');
% Show top 5 regulators for first 3 genes
for g = 1:min(3, length(result.topRegulators))
    fprintf('%s regulators: ', result.topRegulators(g).targetGene);
    fprintf('%s ', result.topRegulators(g).regulatorNames{:});
    fprintf('\n');
end

% Plot GRN
scReservoir_plot_GRN(result.GRN, geneNames, ...
    'threshold', 0.1, ...
    'topN',      40, ...
    'show_net',  false);

%% =====================================================================
%  8. Scalable Mode (for large datasets)
%% =====================================================================
fprintf('\n--- Scalable Mode (sparse reservoir + randomized SVD) ---\n');
fprintf('For datasets >50k cells, use scReservoir_scalable directly:\n\n');
fprintf('  [GRN, topRegs, Hlatent] = scReservoir_scalable(X, pseudotime, geneNames, ...\n');
fprintf('      ''N_res'', 800, ''rankSVD'', 80, ''mode'', ''velocity'');\n\n');

%% =====================================================================
%  9. Save Results
%% =====================================================================
save('results_landscape.mat', 'result', 'geneNames', 'pseudotime', 'fate_labels');
fprintf('Results saved to results_landscape.mat\n');
fprintf('\n=== Pipeline Complete ===\n');
