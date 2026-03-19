%% example02_graph_and_velocity.m
%
% scReservoir: Graph-based reservoir GRN + velocity-mode inference
%
% Demonstrates two advanced methods:
%   A) scGraphResNet - propagates reservoir states across the kNN cell graph
%      (better for branching trajectories)
%   B) Velocity mode GRN - models gene regulation as dx/dt = F(x)
%      (better captures continuous-time regulatory dynamics)

clear; clc; close all;
addpath(genpath('../scReservoir'));

%% =====================================================================
%  1. Load / Simulate Branching Data
%% =====================================================================
rng(1);
nCells   = 800;
nGenes   = 300;

% Simulate a branching trajectory (two lineages)
branch = rand(nCells, 1) > 0.5;     % binary fate assignment
pseudotime = rand(nCells, 1);

% Expression differs by branch after pt = 0.5
X_raw = zeros(nCells, nGenes);
for c = 1:nCells
    pt = pseudotime(c);
    base = randn(1, nGenes);
    if branch(c) == 1
        X_raw(c, :) = base + pt * (1 + randn(1, nGenes) * 0.3);
    else
        X_raw(c, :) = base - pt * (0.5 + randn(1, nGenes) * 0.3);
    end
end
X_raw = max(0, X_raw);    % Non-negative counts

geneNames = arrayfun(@(i) sprintf('Gene%04d', i), 1:nGenes, 'UniformOutput', false);

%% =====================================================================
%  2. Preprocess
%% =====================================================================
[X, geneNames, ~] = scReservoir_preprocess(X_raw, geneNames, ...
    'nHVG', 200, 'normalize', true, 'log_transform', true);
[nCells, nGenes] = size(X);

%% =====================================================================
%  3. Initialize Reservoir
%% =====================================================================
res = scReservoir_init(nGenes, ...
    'N_res',          400, ...
    'rho',            0.9, ...
    'leak_rate',      0.3, ...
    'sparse_density', 0.1, ...
    'lambda',         1e-3, ...
    'seed',           7);

%% =====================================================================
%  4a. Graph Reservoir GRN (scGraphResNet)
%% =====================================================================
fprintf('\n=== Graph Reservoir GRN ===\n');
[GRN_graph, topRegs_graph, H_graph] = scReservoir_graphGRN(res, X, pseudotime, geneNames, ...
    'kNN',    15, ...
    'n_iter', 5, ...
    'k',      5);

fprintf('Top regulators for %s (graph):\n', topRegs_graph(1).targetGene);
for i = 1:5
    fprintf('  %d. %s  (%.3f)\n', i, topRegs_graph(1).regulatorNames{i}, ...
        topRegs_graph(1).scores(i));
end

%% =====================================================================
%  4b. Velocity Mode GRN (scResODE approximation)
%% =====================================================================
fprintf('\n=== Velocity Mode GRN (dx/dt) ===\n');
[GRN_vel, topRegs_vel] = scReservoir_GRN(res, X, pseudotime, geneNames, ...
    'mode',     'velocity', ...
    'k',        5, ...
    'ensemble', 5);

fprintf('Top regulators for %s (velocity):\n', topRegs_vel(1).targetGene);
for i = 1:5
    fprintf('  %d. %s  (%.3f)\n', i, topRegs_vel(1).regulatorNames{i}, ...
        topRegs_vel(1).scores(i));
end

%% =====================================================================
%  5. Compare Graph vs Velocity GRN
%% =====================================================================
figure('Name','Graph vs Velocity GRN','Position',[100 100 1200 500]);

subplot(1, 2, 1);
imagesc(GRN_graph(1:50, 1:50));
colormap(hot); colorbar; title('Graph Reservoir GRN (top 50 genes)');
xlabel('Target'); ylabel('Regulator');

subplot(1, 2, 2);
imagesc(GRN_vel(1:50, 1:50));
colormap(hot); colorbar; title('Velocity GRN (top 50 genes)');
xlabel('Target'); ylabel('Regulator');

sgtitle('scReservoir: GRN Comparison', 'FontWeight','bold');

%% =====================================================================
%  6. Visualize Latent Reservoir Embedding
%% =====================================================================
figure('Name','Reservoir Latent Embedding','Position',[100 100 700 600]);
[~, score] = pca(H_graph, 'NumComponents', 2);
gscatter(score(:,1), score(:,2), branch, [], '.', 15);
xlabel('PC1'); ylabel('PC2');
title('Latent Reservoir States (Graph Reservoir) colored by branch');
legend({'Branch 0','Branch 1'});

%% =====================================================================
%  7. Save
%% =====================================================================
save('results_graph_velocity.mat', 'GRN_graph', 'GRN_vel', ...
     'topRegs_graph', 'topRegs_vel', 'H_graph', 'geneNames');
fprintf('\nResults saved to results_graph_velocity.mat\n');
