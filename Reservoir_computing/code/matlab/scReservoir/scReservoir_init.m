function reservoir = scReservoir_init(nGenes, varargin)
% SCRESERVOIR_INIT  Initialize a random reservoir for single-cell analysis.
%
%   reservoir = scReservoir_init(nGenes)
%   reservoir = scReservoir_init(nGenes, 'Name', Value, ...)
%
%   INPUTS
%     nGenes       Number of input genes (features)
%
%   OPTIONAL NAME-VALUE PAIRS
%     'N_res'        Reservoir size            (default: 500)
%     'rho'          Spectral radius           (default: 0.9)
%     'leak_rate'    Leak/alpha rate [0,1]     (default: 0.3)
%     'input_scale'  Input weight scale        (default: 0.5)
%     'lambda'       Ridge regularization      (default: 1e-3)
%     'sparse_density'  W_res sparsity [0,1]   (default: 0.1; 0 = dense)
%     'seed'         Random seed               (default: 42)
%
%   OUTPUT
%     reservoir      Struct with fields:
%                    .W_res, .W_in, .N_res, .rho, .leak_rate, .lambda
%
%   EXAMPLE
%     res = scReservoir_init(2000, 'N_res', 600, 'rho', 0.9);
%
%   See also: scReservoir_run, scReservoir_GRN, scReservoir_landscape

% --- Parse inputs --------------------------------------------------------
p = inputParser;
addRequired(p,  'nGenes',         @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'N_res',          500,   @isnumeric);
addParameter(p, 'rho',            0.9,   @isnumeric);
addParameter(p, 'leak_rate',      0.3,   @isnumeric);
addParameter(p, 'input_scale',    0.5,   @isnumeric);
addParameter(p, 'lambda',         1e-3,  @isnumeric);
addParameter(p, 'sparse_density', 0.1,   @isnumeric);
addParameter(p, 'seed',           42,    @isnumeric);
parse(p, nGenes, varargin{:});
opt = p.Results;

% --- Reproducibility -----------------------------------------------------
rng(opt.seed);

N_res    = opt.N_res;
nGenes   = opt.nGenes;

% --- Recurrent weight matrix W_res ---------------------------------------
if opt.sparse_density > 0 && opt.sparse_density < 1
    W_res = sprandn(N_res, N_res, opt.sparse_density);
else
    W_res = randn(N_res, N_res) * 0.1;
end

% Normalize spectral radius
eigOpts.tol = 1e-6;
eigOpts.maxit = 300;
try
    lambda_max = max(abs(eigs(W_res, 1, 'LM', eigOpts)));
catch
    lambda_max = max(abs(eig(full(W_res))));
end
if lambda_max > 0
    W_res = W_res * (opt.rho / lambda_max);
end

% --- Input weight matrix W_in --------------------------------------------
W_in = (rand(N_res, nGenes) - 0.5) * 2 * opt.input_scale;

% --- Pack struct ---------------------------------------------------------
reservoir.W_res       = W_res;
reservoir.W_in        = W_in;
reservoir.N_res       = N_res;
reservoir.nGenes      = nGenes;
reservoir.rho         = opt.rho;
reservoir.leak_rate   = opt.leak_rate;
reservoir.lambda      = opt.lambda;
reservoir.input_scale = opt.input_scale;
reservoir.seed        = opt.seed;

fprintf('[scReservoir] Initialized: N_res=%d, rho=%.2f, leak=%.2f, sparse=%.0f%%\n', ...
    N_res, opt.rho, opt.leak_rate, opt.sparse_density*100);
end
