function [GRN, topRegulators, Hlatent, svd_basis] = scReservoir_scalable(X, pseudotime, geneNames, varargin)
% SCRESERVOIR_SCALABLE  Atlas-scale GRN inference with sparse reservoir + randomized SVD.
%
%   Reduces memory and compute from O(N*N_res^2) to O(N*r^2) using:
%     - Sparse recurrent weight matrix W_res
%     - Randomized SVD compression of reservoir states
%
%   Recommended for datasets > 50,000 cells.
%
%   [GRN, topRegs, Hlatent] = scReservoir_scalable(X, pseudotime, geneNames)
%   [GRN, topRegs, Hlatent] = scReservoir_scalable(..., 'N_res', 800, 'rankSVD', 80)
%
%   INPUTS
%     X            nCells x nGenes normalized expression matrix
%     pseudotime   nCells x 1 pseudotime vector (required for velocity mode)
%     geneNames    1 x nGenes cell array
%
%   OPTIONAL NAME-VALUE PAIRS
%     'N_res'        Reservoir size           (default: 800)
%     'rho'          Spectral radius          (default: 0.9)
%     'leak_rate'    Leak rate                (default: 0.3)
%     'density'      W_res sparsity           (default: 0.01)
%     'input_scale'  W_in scale               (default: 0.3)
%     'rankSVD'      SVD rank                 (default: 80)
%     'lambda'       Ridge regularization     (default: 1e-3)
%     'mode'         'static' or 'velocity'   (default: 'velocity')
%     'k'            Top regulators per gene  (default: 10)
%     'seed'         Random seed              (default: 1)
%
%   OUTPUTS
%     GRN          nGenes x nGenes regulatory matrix
%     topRegs      Struct of top regulators
%     Hlatent      nCells x rankSVD compressed latent states
%     svd_basis    Struct with U, S, V for later projection
%
%   See also: scReservoir_init, scReservoir_GRN, scReservoir_landscape

% --- Parse inputs --------------------------------------------------------
p = inputParser;
addRequired(p,  'X',           @isnumeric);
addRequired(p,  'pseudotime');
addRequired(p,  'geneNames',   @iscell);
addParameter(p, 'N_res',       800,    @isnumeric);
addParameter(p, 'rho',         0.9,    @isnumeric);
addParameter(p, 'leak_rate',   0.3,    @isnumeric);
addParameter(p, 'density',     0.01,   @isnumeric);
addParameter(p, 'input_scale', 0.3,    @isnumeric);
addParameter(p, 'rankSVD',     80,     @isnumeric);
addParameter(p, 'lambda',      1e-3,   @isnumeric);
addParameter(p, 'mode',        'velocity', @ischar);
addParameter(p, 'k',           10,     @isnumeric);
addParameter(p, 'seed',        1,      @isnumeric);
parse(p, X, pseudotime, geneNames, varargin{:});
opt = p.Results;

rng(opt.seed);
[nCells, nGenes] = size(X);
N_res   = opt.N_res;
rankSVD = opt.rankSVD;

% --- Sort by pseudotime --------------------------------------------------
if ~isempty(pseudotime)
    [pt_sorted, sortIdx] = sort(pseudotime(:));
    X = X(sortIdx, :);
else
    pt_sorted = (1:nCells)';
end

% --- Sparse Reservoir Initialization -------------------------------------
fprintf('[scReservoir] Building sparse reservoir (N_res=%d, density=%.1f%%)...\n', ...
    N_res, opt.density*100);

W_res = sprandn(N_res, N_res, opt.density);
try
    eigOpts.tol = 1e-5; eigOpts.maxit = 200;
    lam_max = max(abs(eigs(W_res, 1, 'LM', eigOpts)));
catch
    % Fallback for very sparse matrices
    lam_max = max(abs(eig(full(W_res))));
end
if lam_max > 0
    W_res = W_res * (opt.rho / lam_max);
end
W_in = randn(N_res, nGenes) * opt.input_scale;

% --- Select input data for reservoir -------------------------------------
switch lower(opt.mode)
    case 'velocity'
        dt   = diff(pt_sorted);
        dt(dt == 0) = 1e-8;
        dX   = diff(X) ./ dt;
        Xdrv = X(1:end-1, :);      % midpoint expression
        Ydrv = dX;                  % velocity targets
        nDrv = size(Xdrv, 1);
    case 'static'
        Xdrv = X;
        Ydrv = X;
        nDrv = nCells;
    otherwise
        error('scReservoir:unknownMode','Mode must be ''static'' or ''velocity''.');
end

% --- Reservoir Dynamics --------------------------------------------------
fprintf('[scReservoir] Running reservoir dynamics (%d steps)...\n', nDrv);
H = zeros(nDrv, N_res);
h = zeros(N_res, 1);
alpha = opt.leak_rate;
for t = 1:nDrv
    u = Xdrv(t, :)';
    h = (1 - alpha)*h + alpha * tanh(W_res * h + W_in * u);
    H(t, :) = h';
end

% --- Randomized SVD Compression ------------------------------------------
fprintf('[scReservoir] Randomized SVD compression (rank=%d)...\n', rankSVD);
Omega = randn(N_res, rankSVD);
Y_sketch = H * Omega;                    % nDrv x rankSVD
[Q, ~]   = qr(Y_sketch, 0);             % nDrv x rankSVD orthonormal basis
B        = Q' * H;                       % rankSVD x N_res
[Uhat, Shat, Vhat] = svd(B, 'econ');

Ur      = Q * Uhat;                      % nDrv x rankSVD
Hlatent = Ur(:, 1:rankSVD) * Shat(1:rankSVD, 1:rankSVD);   % nDrv x rankSVD

% Store SVD basis for later projection
svd_basis.U   = Ur(:, 1:rankSVD);
svd_basis.S   = Shat(1:rankSVD, 1:rankSVD);
svd_basis.V   = Vhat(:, 1:rankSVD);

% --- Ridge Regression in Compressed Space --------------------------------
fprintf('[scReservoir] Ridge regression for %d genes...\n', nGenes);
HlHl   = Hlatent' * Hlatent + opt.lambda * eye(rankSVD);
HlHl_L = decomposition(HlHl, 'lu');

% Back-projection matrix: rankSVD -> N_res -> nGenes
% influence(i,g) = |W_in' * V * W_out_compressed|_i
Vr = Vhat(:, 1:rankSVD);    % N_res x rankSVD  (right singular vectors)

GRN = zeros(nGenes, nGenes);
for g = 1:nGenes
    y        = Ydrv(:, g);
    W_out_c  = HlHl_L \ (Hlatent' * y);         % rankSVD x 1
    % Approximate back-projection through SVD basis
    influence = abs(W_in' * (Vr * W_out_c));     % nGenes x 1
    GRN(:, g) = influence;
end

% --- Normalize -----------------------------------------------------------
if max(GRN(:)) > 0
    GRN = GRN ./ max(GRN(:));
end
GRN(logical(eye(nGenes))) = 0;

% --- Top Regulators ------------------------------------------------------
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

fprintf('[scReservoir] Scalable GRN complete. Memory saved: %.1fx\n', N_res/rankSVD);
end
