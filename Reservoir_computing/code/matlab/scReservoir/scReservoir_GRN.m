function [GRN, topRegulators] = scReservoir_GRN(reservoir, X, pseudotime, geneNames, varargin)
% SCRESERVOIR_GRN  Infer gene regulatory network from single-cell data.
%
%   [GRN, topRegs] = scReservoir_GRN(reservoir, X, pseudotime, geneNames)
%   [GRN, topRegs] = scReservoir_GRN(..., 'mode', 'causal', 'k', 10)
%
%   INPUTS
%     reservoir    Struct from scReservoir_init
%     X            nCells x nGenes expression matrix
%     pseudotime   nCells x 1 pseudotime vector (use [] to skip ordering)
%     geneNames    1 x nGenes cell array of gene name strings
%
%   OPTIONAL NAME-VALUE PAIRS
%     'mode'       Inference mode:
%                  'static'   - regress expression on reservoir states (default)
%                  'causal'   - lag reservoir states for directional causality
%                  'velocity' - regress gene velocity (dx/dt) on reservoir states
%     'k'          Number of top regulators to return per gene (default: 10)
%     'threshold'  GRN edge threshold [0,1] (default: 0.1)
%     'normalize'  Normalize GRN to [0,1]   (default: true)
%     'ensemble'   Number of reservoir ensembles (default: 1)
%
%   OUTPUTS
%     GRN              nGenes x nGenes normalized influence matrix
%                      GRN(i,j) = influence of gene i on gene j
%     topRegulators    Struct array (1 x nGenes) with fields:
%                      .targetGene, .regulatorNames, .regulatorIdx, .scores
%
%   EXAMPLE
%     [GRN, topRegs] = scReservoir_GRN(res, X, pt, genes, 'mode','causal','k',5);
%     disp(topRegs(10))
%
%   See also: scReservoir_init, scReservoir_run, scReservoir_graphGRN

% --- Parse inputs --------------------------------------------------------
p = inputParser;
addRequired(p, 'reservoir');
addRequired(p, 'X',           @isnumeric);
addRequired(p, 'pseudotime');
addRequired(p, 'geneNames',   @iscell);
addParameter(p,'mode',        'static',  @ischar);
addParameter(p,'k',           10,        @isnumeric);
addParameter(p,'threshold',   0.1,       @isnumeric);
addParameter(p,'normalize',   true,      @islogical);
addParameter(p,'ensemble',    1,         @isnumeric);
parse(p, reservoir, X, pseudotime, geneNames, varargin{:});
opt = p.Results;

[nCells, nGenes] = size(X);
lambda = reservoir.lambda;

% --- Sort cells by pseudotime --------------------------------------------
if ~isempty(pseudotime)
    [~, sortIdx] = sort(pseudotime(:));
    X = X(sortIdx, :);
end

% --- Ensemble averaging --------------------------------------------------
GRN_all = zeros(nGenes, nGenes, opt.ensemble);

for e = 1:opt.ensemble
    if opt.ensemble > 1
        % Re-initialize reservoir with different seed for ensemble
        res_e = scReservoir_init(nGenes, ...
            'N_res',          reservoir.N_res, ...
            'rho',            reservoir.rho, ...
            'leak_rate',      reservoir.leak_rate, ...
            'input_scale',    reservoir.input_scale, ...
            'lambda',         reservoir.lambda, ...
            'seed',           reservoir.seed + e - 1);
    else
        res_e = reservoir;
    end

    W_res = res_e.W_res;
    W_in  = res_e.W_in;
    N_res = res_e.N_res;
    alpha = res_e.leak_rate;

    % --- Reservoir dynamics ----------------------------------------------
    H = zeros(nCells, N_res);
    h = zeros(N_res, 1);
    for t = 1:nCells
        u = X(t, :)';
        h = (1 - alpha)*h + alpha * tanh(W_res*h + W_in*u);
        H(t, :) = h';
    end

    % --- Mode-specific regression input preparation ----------------------
    switch lower(opt.mode)
        case 'static'
            H_reg = H;
            Y_reg = X;

        case 'causal'
            % Lag reservoir by one step to enforce temporal causality
            H_reg = [zeros(1, N_res); H(1:end-1, :)];
            Y_reg = X;

        case 'velocity'
            % Regress gene velocity dx/dt instead of expression
            if isempty(pseudotime)
                error('scReservoir:noPseudotime', ...
                    'Velocity mode requires pseudotime input.');
            end
            pt_sorted = sort(pseudotime(:));
            dt = diff(pt_sorted);                      % (nCells-1) x 1
            dt(dt == 0) = 1e-8;                        % avoid division by zero
            dX = diff(X) ./ dt;                        % (nCells-1) x nGenes
            X_mid = X(1:end-1, :);

            % Re-compute reservoir on midpoints
            H_mid = zeros(nCells-1, N_res);
            h = zeros(N_res, 1);
            for t = 1:nCells-1
                u = X_mid(t, :)';
                h = (1 - alpha)*h + alpha * tanh(W_res*h + W_in*u);
                H_mid(t, :) = h';
            end
            H_reg = H_mid;
            Y_reg = dX;

        otherwise
            error('scReservoir:unknownMode', ...
                'Unknown mode "%s". Use ''static'', ''causal'', or ''velocity''.', opt.mode);
    end

    % --- Infer GRN via ridge regression ----------------------------------
    HtH    = H_reg' * H_reg + lambda * eye(N_res);
    HtH_LU = decomposition(HtH, 'lu');   % Reuse factorization

    GRN_e = zeros(nGenes, nGenes);
    for g = 1:nGenes
        y      = Y_reg(:, g);
        W_out  = HtH_LU \ (H_reg' * y);          % Ridge readout (N_res x 1)
        influence = abs(W_in' * W_out);            % Back-project to genes (nGenes x 1)
        GRN_e(:, g) = influence;
    end

    GRN_all(:, :, e) = GRN_e;
end

% --- Average ensemble ----------------------------------------------------
GRN = mean(GRN_all, 3);

% --- Normalize -----------------------------------------------------------
if opt.normalize && max(GRN(:)) > 0
    GRN = GRN ./ max(GRN(:));
end

% Remove self-regulation
GRN(logical(eye(nGenes))) = 0;

% --- Extract top regulators ----------------------------------------------
topRegulators = struct('targetGene',     cell(1, nGenes), ...
                       'regulatorNames', cell(1, nGenes), ...
                       'regulatorIdx',   cell(1, nGenes), ...
                       'scores',         cell(1, nGenes));
k = min(opt.k, nGenes - 1);

for g = 1:nGenes
    col = GRN(:, g);
    col(g) = -inf;   % exclude self
    [vals, idx] = sort(col, 'descend');
    topRegulators(g).targetGene     = geneNames{g};
    topRegulators(g).regulatorIdx   = idx(1:k);
    topRegulators(g).regulatorNames = geneNames(idx(1:k));
    topRegulators(g).scores         = vals(1:k);
end

fprintf('[scReservoir] GRN inferred: mode=%s, ensemble=%d, nGenes=%d\n', ...
    opt.mode, opt.ensemble, nGenes);
end
