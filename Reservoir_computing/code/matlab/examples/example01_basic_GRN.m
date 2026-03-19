%% example01_basic_GRN.m
%
% scReservoir: Basic GRN inference from scRNA-seq data
%
% This example demonstrates the core scReservoir workflow:
%   1. Preprocess a normalized expression matrix
%   2. Initialize a random reservoir
%   3. Infer a gene regulatory network in static and causal modes
%   4. Display top regulators and visualize the GRN
%
% DATA FORMAT:
%   X         - nCells x nGenes raw or normalized expression matrix
%   geneNames - 1 x nGenes cell array of gene name strings
%   pseudotime- nCells x 1 numeric vector (e.g., from Monocle or DPT)
%
% Replace the synthetic data block below with your own data.

clear; clc; close all;
addpath(genpath('../scReservoir'));

%% =====================================================================
%  1. Load / Simulate Data
%% =====================================================================

% --- Synthetic data (replace with your own) ---
rng(0);
nCells  = 500;
nGenes  = 200;

% Simulated pseudotime
pseudotime = linspace(0, 1, nCells)' + 0.05*randn(nCells, 1);

% Simulated gene expression with some regulatory structure
X_raw = max(0, randn(nCells, nGenes) + repmat(pseudotime, 1, nGenes) .* randn(1, nGenes));

% Gene names (replace with your actual gene names)
geneNames = arrayfun(@(i) sprintf('Gene%04d', i), 1:nGenes, 'UniformOutput', false);

fprintf('Data: %d cells x %d genes\n', nCells, nGenes);

%% =====================================================================
%  2. Preprocess
%% =====================================================================
[X, geneNames, ~] = scReservoir_preprocess(X_raw, geneNames, ...
    'nHVG',    150, ...       % Select 150 highly variable genes
    'normalize', true, ...
    'log_transform', true);

[nCells, nGenes] = size(X);
fprintf('After HVG selection: %d cells x %d genes\n', nCells, nGenes);

%% =====================================================================
%  3. Initialize Reservoir
%% =====================================================================
res = scReservoir_init(nGenes, ...
    'N_res',          300, ...     % Reservoir size (increase for real data)
    'rho',            0.9, ...     % Spectral radius
    'leak_rate',      0.3, ...     % Temporal smoothing
    'sparse_density', 0.1, ...     % 10% sparsity
    'lambda',         1e-3, ...    % Ridge regularization
    'seed',           42);

%% =====================================================================
%  4a. Static GRN (no causality)
%% =====================================================================
fprintf('\n--- Static GRN ---\n');
[GRN_static, topRegs_static] = scReservoir_GRN(res, X, pseudotime, geneNames, ...
    'mode', 'static', ...
    'k',    5);

% Display top regulators for gene 1
fprintf('\nTop regulators for %s (static):\n', topRegs_static(1).targetGene);
for i = 1:5
    fprintf('  %d. %s  (score=%.3f)\n', i, ...
        topRegs_static(1).regulatorNames{i}, ...
        topRegs_static(1).scores(i));
end

%% =====================================================================
%  4b. Causal GRN (pseudotime-ordered, directional)
%% =====================================================================
fprintf('\n--- Causal GRN (pseudotime-ordered) ---\n');
[GRN_causal, topRegs_causal] = scReservoir_GRN(res, X, pseudotime, geneNames, ...
    'mode', 'causal', ...
    'k',    5, ...
    'ensemble', 3);   % Average over 3 reservoir initializations

% Display top regulators for gene 1
fprintf('\nTop regulators for %s (causal):\n', topRegs_causal(1).targetGene);
for i = 1:5
    fprintf('  %d. %s  (score=%.3f)\n', i, ...
        topRegs_causal(1).regulatorNames{i}, ...
        topRegs_causal(1).scores(i));
end

%% =====================================================================
%  5. Visualize GRN
%% =====================================================================
scReservoir_plot_GRN(GRN_causal, geneNames, ...
    'threshold', 0.15, ...
    'topN',      30, ...
    'show_net',  true);

%% =====================================================================
%  6. Save results
%% =====================================================================
save('results_basic_GRN.mat', 'GRN_static', 'GRN_causal', ...
     'topRegs_static', 'topRegs_causal', 'geneNames', 'res');
fprintf('\nResults saved to results_basic_GRN.mat\n');
