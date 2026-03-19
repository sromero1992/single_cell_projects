function result = scReservoir_landscape(X, pseudotime, geneNames, varargin)
% SCRESERVOIR_LANDSCAPE  Reconstruct Waddington landscape and infer cell fate attractors.
%
%   Full pipeline: sparse reservoir -> randomized SVD -> latent dynamics ->
%   attractor detection -> energy landscape -> fate probabilities -> GRN.
%
%   result = scReservoir_landscape(X, pseudotime, geneNames)
%   result = scReservoir_landscape(..., 'nAttractors', 5, 'N_res', 600)
%
%   INPUTS
%     X            nCells x nGenes normalized expression matrix
%     pseudotime   nCells x 1 pseudotime vector
%     geneNames    1 x nGenes cell array of gene names
%
%   OPTIONAL NAME-VALUE PAIRS
%     'N_res'         Reservoir size              (default: 600)
%     'rho'           Spectral radius             (default: 0.9)
%     'leak_rate'     Leak rate                   (default: 0.3)
%     'density'       W_res sparsity              (default: 0.01)
%     'input_scale'   W_in scale                  (default: 0.3)
%     'rankSVD'       SVD rank for compression    (default: 50)
%     'lambda'        Ridge regularization        (default: 1e-3)
%     'nAttractors'   Number of attractor states  (default: 5)
%     'attractor_pct' Velocity percentile cutoff  (default: 5)
%     'k'             Top regulators per gene     (default: 10)
%     'seed'          Random seed                 (default: 1)
%     'plot'          Generate diagnostic plots   (default: true)
%
%   OUTPUT  result struct with fields:
%     .Hlatent          nCells x rankSVD latent reservoir states
%     .A                rankSVD x rankSVD latent dynamics matrix
%     .velocity         nCells x rankSVD velocity field
%     .speed            nCells x 1 velocity magnitudes
%     .attractorCells   indices of low-velocity (attractor) cells
%     .attractorLabels  cluster labels for attractor cells
%     .attractorCenters nAttractors x rankSVD cluster centroids
%     .attractorGenes   nAttractors x nGenes mean expression per attractor
%     .fateProbabilities nCells x nAttractors soft assignment probabilities
%     .GRN              nGenes x nGenes regulatory matrix
%     .topRegulators    struct array of top regulators per gene
%     .energy           nCells x 1 approximate Waddington potential
%     .pca_scores       nCells x 2 PCA projection of Hlatent
%
%   See also: scReservoir_init, scReservoir_scalable, scReservoir_plot

% --- Parse inputs --------------------------------------------------------
p = inputParser;
addRequired(p,  'X',             @isnumeric);
addRequired(p,  'pseudotime',    @isnumeric);
addRequired(p,  'geneNames',     @iscell);
addParameter(p, 'N_res',         600,   @isnumeric);
addParameter(p, 'rho',           0.9,   @isnumeric);
addParameter(p, 'leak_rate',     0.3,   @isnumeric);
addParameter(p, 'density',       0.01,  @isnumeric);
addParameter(p, 'input_scale',   0.3,   @isnumeric);
addParameter(p, 'rankSVD',       50,    @isnumeric);
addParameter(p, 'lambda',        1e-3,  @isnumeric);
addParameter(p, 'nAttractors',   5,     @isnumeric);
addParameter(p, 'attractor_pct', 5,     @isnumeric);
addParameter(p, 'k',             10,    @isnumeric);
addParameter(p, 'seed',          1,     @isnumeric);
addParameter(p, 'plot',          true,  @islogical);
parse(p, X, pseudotime, geneNames, varargin{:});
opt = p.Results;

rng(opt.seed);
[nCells, nGenes] = size(X);
N_res   = opt.N_res;
rankSVD = opt.rankSVD;

%% =====================================================================
%  STEP 1: Sort cells by pseudotime
%% =====================================================================
[pt_sorted, sortIdx] = sort(pseudotime(:));
X = X(sortIdx, :);

%% =====================================================================
%  STEP 2: Sparse Reservoir Initialization
%% =====================================================================
fprintf('[scReservoir] Initializing sparse reservoir...\n');
W_res = sprandn(N_res, N_res, opt.density);
try
    eigOpts.tol = 1e-5; eigOpts.maxit = 200;
    lam_max = max(abs(eigs(W_res, 1, 'LM', eigOpts)));
catch
    lam_max = max(abs(eig(full(W_res))));
end
if lam_max > 0, W_res = W_res * (opt.rho / lam_max); end
W_in = randn(N_res, nGenes) * opt.input_scale;

%% =====================================================================
%  STEP 3: Reservoir Dynamics
%% =====================================================================
fprintf('[scReservoir] Running reservoir dynamics...\n');
H = zeros(nCells, N_res);
h = zeros(N_res, 1);
alpha = opt.leak_rate;
for t = 1:nCells
    u = X(t, :)';
    h = (1 - alpha)*h + alpha * tanh(W_res*h + W_in*u);
    H(t, :) = h';
end

%% =====================================================================
%  STEP 4: Randomized SVD Compression
%% =====================================================================
fprintf('[scReservoir] Randomized SVD (rank=%d)...\n', rankSVD);
Omega    = randn(N_res, rankSVD);
Y_sketch = H * Omega;
[Q, ~]   = qr(Y_sketch, 0);
B        = Q' * H;
[Uhat, Shat, Vhat] = svd(B, 'econ');
Ur      = Q * Uhat;
Hlatent = Ur(:, 1:rankSVD) * Shat(1:rankSVD, 1:rankSVD);

%% =====================================================================
%  STEP 5: Latent Dynamics Estimation  (dh/dt = A*h)
%% =====================================================================
fprintf('[scReservoir] Estimating latent dynamics...\n');
dt  = diff(pt_sorted);
dt(dt == 0) = 1e-8;
dH  = diff(Hlatent) ./ dt;
Hmid = Hlatent(1:end-1, :);

A = (Hmid' * Hmid + opt.lambda * eye(rankSVD)) \ (Hmid' * dH);  % rankSVD x rankSVD

%% =====================================================================
%  STEP 6: Velocity Field and Attractor Detection
%% =====================================================================
fprintf('[scReservoir] Computing velocity field...\n');
velocity = Hlatent * A';                         % nCells x rankSVD
speed    = sqrt(sum(velocity .^ 2, 2));          % nCells x 1

threshold      = prctile(speed, opt.attractor_pct);
attractorCells = find(speed < threshold);

fprintf('[scReservoir] Clustering %d attractor cells into %d fates...\n', ...
    numel(attractorCells), opt.nAttractors);
[attractorLabels, attractorCenters] = kmeans( ...
    Hlatent(attractorCells, :), opt.nAttractors, ...
    'Replicates', 5, 'MaxIter', 500);

%% =====================================================================
%  STEP 7: Gene Programs Per Attractor
%% =====================================================================
attractorGenes = zeros(opt.nAttractors, nGenes);
for i = 1:opt.nAttractors
    cellsInCluster = attractorCells(attractorLabels == i);
    if ~isempty(cellsInCluster)
        attractorGenes(i, :) = mean(X(cellsInCluster, :), 1);
    end
end

%% =====================================================================
%  STEP 8: Fate Probabilities (softmax over distance to attractor centers)
%% =====================================================================
fprintf('[scReservoir] Computing fate probabilities...\n');
distToAttractors = pdist2(Hlatent, attractorCenters);   % nCells x nAttractors
prob = exp(-distToAttractors);
prob = prob ./ sum(prob, 2);                             % Row normalize

%% =====================================================================
%  STEP 9: GRN from Latent Dynamics  (Gamma = |W_in' * A * W_in|)
%% =====================================================================
fprintf('[scReservoir] Computing GRN from latent dynamics...\n');
% Map: latent dynamics A (rankSVD x rankSVD) back to gene space via W_in
% Full back-projection: W_in' * V_r * A * V_r' * W_in (approx)
Vr  = Vhat(:, 1:rankSVD);           % N_res x rankSVD
GRN = abs(W_in' * (Vr * A * Vr') * W_in);   % nGenes x nGenes
if max(GRN(:)) > 0
    GRN = GRN ./ max(GRN(:));
end
GRN(logical(eye(nGenes))) = 0;

% Top regulators
k = min(opt.k, nGenes - 1);
topRegulators = struct('targetGene',     cell(1, nGenes), ...
                       'regulatorNames', cell(1, nGenes), ...
                       'regulatorIdx',   cell(1, nGenes), ...
                       'scores',         cell(1, nGenes));
for g = 1:nGenes
    col = GRN(:, g); col(g) = -inf;
    [vals, idx] = sort(col, 'descend');
    topRegulators(g).targetGene     = geneNames{g};
    topRegulators(g).regulatorIdx   = idx(1:k);
    topRegulators(g).regulatorNames = geneNames(idx(1:k));
    topRegulators(g).scores         = vals(1:k);
end

%% =====================================================================
%  STEP 10: Energy Landscape (Waddington potential)
%% =====================================================================
fprintf('[scReservoir] Reconstructing energy landscape...\n');
[~, pca_scores] = pca(Hlatent, 'NumComponents', 2);
pc1 = pca_scores(:, 1);
pc2 = pca_scores(:, 2);

% KDE-based density estimate for potential energy
[density, ~] = ksdensity([pc1, pc2], [pc1, pc2]);
density = max(density, 1e-10);    % avoid log(0)
energy  = -log(density);
energy  = energy - min(energy);   % shift to start at 0

%% =====================================================================
%  STEP 11: Visualization
%% =====================================================================
if opt.plot
    scReservoir_plot_landscape(pc1, pc2, energy, prob, pseudotime(sortIdx), attractorCells);
end

%% =====================================================================
%  STEP 12: Pack output struct
%% =====================================================================
result.Hlatent           = Hlatent;
result.A                 = A;
result.velocity          = velocity;
result.speed             = speed;
result.attractorCells    = attractorCells;
result.attractorLabels   = attractorLabels;
result.attractorCenters  = attractorCenters;
result.attractorGenes    = attractorGenes;
result.fateProbabilities = prob;
result.GRN               = GRN;
result.topRegulators     = topRegulators;
result.energy            = energy;
result.pca_scores        = pca_scores;
result.W_in              = W_in;
result.W_res             = W_res;
result.A_latent          = A;

fprintf('[scReservoir] Pipeline complete.\n');
fprintf('  Attractors found: %d  |  Genes: %d  |  Latent rank: %d\n', ...
    opt.nAttractors, nGenes, rankSVD);
end
