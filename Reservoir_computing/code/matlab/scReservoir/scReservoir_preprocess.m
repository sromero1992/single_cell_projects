function [X_out, geneNames_out, hvgIdx] = scReservoir_preprocess(X, geneNames, varargin)
% SCRESERVOIR_PREPROCESS  Normalize and select highly variable genes.
%
%   [X_out, genes_out, hvgIdx] = scReservoir_preprocess(X, geneNames)
%   [X_out, genes_out, hvgIdx] = scReservoir_preprocess(X, geneNames, 'nHVG', 2000)
%
%   INPUTS
%     X          nCells x nGenes raw count matrix
%     geneNames  1 x nGenes cell array of gene names
%
%   OPTIONAL NAME-VALUE PAIRS
%     'nHVG'          Number of highly variable genes to select (default: 2000)
%                     Set to Inf to keep all genes.
%     'normalize'     Library-size normalize to target_sum (default: true)
%     'target_sum'    Target counts per cell (default: 1e4)
%     'log_transform' Apply log1p transform (default: true)
%     'scale'         Z-score scale genes after log (default: false)
%     'min_cells'     Minimum cells expressing a gene to keep it (default: 3)
%     'min_counts'    Minimum total counts per gene (default: 1)
%
%   OUTPUTS
%     X_out       nCells x nHVG processed matrix
%     genes_out   1 x nHVG filtered gene names
%     hvgIdx      Indices into original geneNames for selected HVGs

p = inputParser;
addRequired(p,  'X',            @isnumeric);
addRequired(p,  'geneNames',    @iscell);
addParameter(p, 'nHVG',         2000,  @isnumeric);
addParameter(p, 'normalize',    true,  @islogical);
addParameter(p, 'target_sum',   1e4,   @isnumeric);
addParameter(p, 'log_transform',true,  @islogical);
addParameter(p, 'scale',        false, @islogical);
addParameter(p, 'min_cells',    3,     @isnumeric);
addParameter(p, 'min_counts',   1,     @isnumeric);
parse(p, X, geneNames, varargin{:});
opt = p.Results;

[nCells, nGenes] = size(X);
fprintf('[scReservoir] Preprocessing: %d cells x %d genes\n', nCells, nGenes);

% --- Filter low-expressed genes -----------------------------------------
n_cells_expressed = sum(X > 0, 1);      % 1 x nGenes
total_counts      = sum(X, 1);          % 1 x nGenes
geneKeep = (n_cells_expressed >= opt.min_cells) & (total_counts >= opt.min_counts);
X = X(:, geneKeep);
geneNames = geneNames(geneKeep);
fprintf('[scReservoir]   After gene filter: %d genes remain\n', sum(geneKeep));

% --- Library-size normalization -----------------------------------------
if opt.normalize
    lib_size = sum(X, 2);                    % nCells x 1
    lib_size(lib_size == 0) = 1;
    X = bsxfun(@rdivide, X, lib_size) * opt.target_sum;
end

% --- Log1p transform -----------------------------------------------------
if opt.log_transform
    X = log1p(X);
end

% --- Highly Variable Gene selection -------------------------------------
nHVG = min(opt.nHVG, size(X, 2));

if isinf(opt.nHVG) || nHVG >= size(X, 2)
    hvgIdx    = 1:size(X, 2);
    X_out     = X;
    geneNames_out = geneNames;
else
    % Compute mean and dispersion (variance/mean) per gene
    gene_mean = mean(X, 1);
    gene_var  = var(X, 0, 1);
    dispersion = gene_var ./ (gene_mean + 1e-8);

    % Bin-normalize dispersion by mean (Seurat-style)
    nBins = 20;
    [~, ~, binIdx] = histcounts(gene_mean, nBins);
    binIdx(binIdx == 0) = 1;
    disp_norm = zeros(1, size(X, 2));
    for b = 1:nBins
        inBin = (binIdx == b);
        if sum(inBin) > 1
            mu_d  = mean(dispersion(inBin));
            std_d = std(dispersion(inBin));
            if std_d > 0
                disp_norm(inBin) = (dispersion(inBin) - mu_d) / std_d;
            end
        end
    end

    [~, sortedDisp] = sort(disp_norm, 'descend');
    hvgLocal  = sortedDisp(1:nHVG);
    hvgIdx    = find(geneKeep);
    hvgIdx    = hvgIdx(hvgLocal);

    X_out     = X(:, hvgLocal);
    geneNames_out = geneNames(hvgLocal);
    fprintf('[scReservoir]   Selected %d HVGs\n', nHVG);
end

% --- Z-score scaling (optional) -----------------------------------------
if opt.scale
    X_out = (X_out - mean(X_out, 1)) ./ (std(X_out, 0, 1) + 1e-8);
    X_out = max(-10, min(10, X_out));   % clip extreme values
end
end
