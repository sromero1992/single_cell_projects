function H = scReservoir_run(reservoir, X, varargin)
% SCRESERVOIR_RUN  Drive reservoir with gene expression input, return state matrix.
%
%   H = scReservoir_run(reservoir, X)
%   H = scReservoir_run(reservoir, X, 'washout', 50)
%
%   INPUTS
%     reservoir    Struct from scReservoir_init
%     X            nCells x nGenes expression matrix (normalized, log-transformed)
%
%   OPTIONAL NAME-VALUE PAIRS
%     'washout'    Number of initial states to discard (default: 0)
%     'h0'         Initial reservoir state N_res x 1   (default: zeros)
%
%   OUTPUT
%     H            nCells x N_res reservoir state matrix
%
%   The leaky integrator update rule is:
%     h(t) = (1-alpha)*h(t-1) + alpha*tanh(W_res*h(t-1) + W_in*x(t))
%
%   See also: scReservoir_init, scReservoir_GRN

p = inputParser;
addRequired(p,  'reservoir');
addRequired(p,  'X',          @isnumeric);
addParameter(p, 'washout',    0,    @isnumeric);
addParameter(p, 'h0',         [],   @isnumeric);
parse(p, reservoir, X, varargin{:});
opt = p.Results;

[nCells, nGenes] = size(X);
N_res  = reservoir.N_res;
alpha  = reservoir.leak_rate;
W_res  = reservoir.W_res;
W_in   = reservoir.W_in;

% Validate dimensions
if nGenes ~= reservoir.nGenes
    error('scReservoir:dimMismatch', ...
        'X has %d genes but reservoir expects %d.', nGenes, reservoir.nGenes);
end

% Initial state
if isempty(opt.h0)
    h = zeros(N_res, 1);
else
    h = opt.h0(:);
end

H = zeros(nCells, N_res);

for t = 1:nCells
    u = X(t, :)';
    h = (1 - alpha) * h + alpha * tanh(W_res * h + W_in * u);
    H(t, :) = h';
end

% Discard washout steps
if opt.washout > 0 && opt.washout < nCells
    H = H(opt.washout+1:end, :);
end
end
