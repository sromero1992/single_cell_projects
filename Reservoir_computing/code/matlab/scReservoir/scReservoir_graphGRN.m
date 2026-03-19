function [GRN, topRegulators, H] = scReservoir_graphGRN(reservoir, X, pseudotime, geneNames, varargin)
% SCRESERVOIR_GRAPHGRN  Graph-based reservoir GRN inference (scGraphResNet).
%
%   Propagates reservoir states across the k-nearest-neighbor cell graph,
%   capturing branching trajectories and nonlinear manifold structure.
%   Better suited than sequential reservoir for non-linear developmental data.
%
%   [GRN, topRegs, H] = scReservoir_graphGRN(reservoir, X, pseudotime, geneNames)
%   [GRN, topRegs, H] = scReservoir_graphGRN(..., 'kNN', 15, 'n_iter', 5)
%
%   INPUTS
%     reservoir    Struct from scReservoir_init
%     X            nCells x nGenes expression matrix
%     pseudotime   nCells x 1 vector (used only for sorting; can be [])
%     geneNames    1 x nGenes cell array
%
%   OPTIONAL NAME-VALUE PAIRS
%     'kNN'        Number of nearest neighbors for cell graph (default: 10)
%     'n_iter'     Graph reservoir propagation iterations   (default: 5)
%     'k'          Top regulators per gene                 (default: 10)
%     'normalize'  Normalize GRN to [0,1]                 (default: true)
%     'dist_metric' Distance for kNN: 'euclidean','cosine' (default: 'euclidean')
%
%   OUTPUTS
%     GRN              nGenes x nGenes regulatory matrix
%     topRegulators    Struct array of top regulators per gene
%     H                nCells x N_res final reservoir states
%
%   See also: scReservoir_init, scReservoir_GRN, scReservoir_landscape

% --- Parse inputs --------------------------------------------------------
p = inputParser;
addRequired(p, 'reservoir');
addRequired(p, 'X',            @isnumeric);
addRequired(p, 'pseudotime');
addRequired(p, 'geneNames',    @iscell);
addParameter(p,'kNN',          10,           @isnumeric);
addParameter(p,'n_iter',       5,            @isnumeric);
addParameter(p,'k',            10,           @isnumeric);
addParameter(p,'normalize',    true,         @islogical);
addParameter(p,'dist_metric',  'euclidean',  @ischar);
parse(p, reservoir, X, pseudotime, geneNames, varargin{:});
opt = p.Results;

[nCells, nGenes] = size(X);
N_res  = reservoir.N_res;
alpha  = reservoir.leak_rate;
lambda = reservoir.lambda;
W_res  = reservoir.W_res;
W_in   = reservoir.W_in;

% --- Sort by pseudotime --------------------------------------------------
if ~isempty(pseudotime)
    [~, sortIdx] = sort(pseudotime(:));
    X = X(sortIdx, :);
end

% --- Build kNN cell graph ------------------------------------------------
fprintf('[scReservoir] Building kNN graph (k=%d)...\n', opt.kNN);
D = pdist2(X, X, opt.dist_metric);
[~, neighbors] = sort(D, 2);     % nCells x nCells, sorted by distance

% Build adjacency (exclude self = column 1)
A = zeros(nCells, nCells);
for i = 1:nCells
    nbrs = neighbors(i, 2:opt.kNN+1);
    A(i, nbrs) = 1;
end
A = max(A, A');                   % Symmetrize

% Row-normalize (random walk normalization)
deg = sum(A, 2);
deg(deg == 0) = 1;
G_norm = bsxfun(@rdivide, A, deg);   % nCells x nCells row-stochastic

% --- Graph Reservoir Propagation -----------------------------------------
fprintf('[scReservoir] Running graph reservoir (%d iterations)...\n', opt.n_iter);
H = zeros(nCells, N_res);

for iter = 1:opt.n_iter
    H_new = zeros(nCells, N_res);
    for i = 1:nCells
        % Aggregate neighbor states (row of G_norm gives weights)
        neighbor_state = (G_norm(i, :) * H)';   % N_res x 1 weighted aggregate
        u = X(i, :)';
        h_new = tanh(W_res * neighbor_state + W_in * u);
        H_new(i, :) = (1 - alpha) * H(i, :) + alpha * h_new';
    end
    H = H_new;
end

% --- GRN Inference -------------------------------------------------------
fprintf('[scReservoir] Inferring GRN...\n');
HtH    = H' * H + lambda * eye(N_res);
HtH_LU = decomposition(HtH, 'lu');

GRN = zeros(nGenes, nGenes);
for g = 1:nGenes
    y         = X(:, g);
    W_out     = HtH_LU \ (H' * y);
    influence = abs(W_in' * W_out);
    GRN(:, g) = influence;
end

% --- Normalize -----------------------------------------------------------
if opt.normalize && max(GRN(:)) > 0
    GRN = GRN ./ max(GRN(:));
end
GRN(logical(eye(nGenes))) = 0;

% --- Top regulators ------------------------------------------------------
k = min(opt.k, nGenes - 1);
topRegulators = struct('targetGene',     cell(1, nGenes), ...
                       'regulatorNames', cell(1, nGenes), ...
                       'regulatorIdx',   cell(1, nGenes), ...
                       'scores',         cell(1, nGenes));
for g = 1:nGenes
    col = GRN(:, g);
    col(g) = -inf;
    [vals, idx] = sort(col, 'descend');
    topRegulators(g).targetGene     = geneNames{g};
    topRegulators(g).regulatorIdx   = idx(1:k);
    topRegulators(g).regulatorNames = geneNames(idx(1:k));
    topRegulators(g).scores         = vals(1:k);
end

fprintf('[scReservoir] Graph GRN complete.\n');
end
