function [Q, Q1, P] = assemble_qubo_matrix(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K, Xnet_target, penalty_scale, auto_penalty)
% ASSEMBLE_QUBO_MATRIX  Construct the QUBO matrix with cardinality penalty.
%
%   [Q, Q1, P] = ASSEMBLE_QUBO_MATRIX(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K)
%   [Q, Q1, P] = ASSEMBLE_QUBO_MATRIX(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K, Xnet_target)
%   [Q, Q1, P] = ASSEMBLE_QUBO_MATRIX(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K, Xnet_target, penalty_scale)
%   [Q, Q1, P] = ASSEMBLE_QUBO_MATRIX(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K, Xnet_target, penalty_scale, auto_penalty)
%
% DESCRIPTION:
%   Assembles the base QUBO matrix
%
%       Q1 = Xdiff + diag(Vdiff) - (Xmnn_A + Xmnn_B)  [+ Xnet_target]
%
%   then adds the cardinality penalty to enforce exactly K gene selections:
%
%       off-diagonal: Q(i,j) = Q1(i,j) + P    (i ~= j)
%       diagonal:     Q(i,i) = Q1(i,i) + P*(1 - 2*K)
%
%   Penalty coefficient P is derived automatically from the spectral norm
%   of Q1 (auto_penalty = true, default):
%
%       P = penalty_scale * norm(Q1, 2)     [spectral norm, = sigma_max(Q1)]
%
%   The spectral norm equals the largest absolute eigenvalue for symmetric
%   Q1 and provides a principled lower bound: any solution violating
%   |z| = K incurs a penalty energy exceeding the maximum biological energy
%   gain, which is bounded by the spectral radius of Q1.  This is tighter
%   than the legacy elementwise heuristic P = penalty_scale * max(|Q1|),
%   which ignores matrix structure and over-penalises on the low-rank
%   co-expression matrices typical of scRNA-seq data.
%
%   Set auto_penalty = false to revert to the legacy heuristic.
%
% INPUTS (required):
%   Xdiff         - Differential co-expression matrix, G×G.
%                   X_diff = S_B - S_A (condition B minus condition A).
%                   Negative off-diagonal entries mean co-expression
%                   increased in condition A (test).
%   Vdiff         - Gene-wise cell-state difference vector, G×1.
%                   V_diff(g) = v_g^(B) - v_g^(A), where v_g is the
%                   projection of gene g's normalised expression onto the
%                   per-cell biological state scalar of that condition.
%                   The scalar can be anything continuous: UCell pathway
%                   activity score, pseudotime, cell potency, differentiation
%                   rank, etc.  Pass zeros(G,1) to omit this term.
%   Xmnn_A        - MNN adjacency matrix for condition A (test),  G×G sparse.
%   Xmnn_B        - MNN adjacency matrix for condition B (reference), G×G sparse.
%   K             - Target subnetwork size (positive integer).
%
% INPUTS (optional):
%   Xnet_target   - Pathway membership prior matrix, G×G.
%                   Off-diagonal entry (i,j) = -1 when both genes are core
%                   pathway members, 0 otherwise.
%                   IMPORTANT: when the analysis is already restricted to
%                   pathway genes (standard use), this matrix is uniformly
%                   -1 off-diagonal and adds only a constant energy shift
%                   that does NOT affect which solution is optimal.
%                   Provide it only when extra candidate genes (TF targets,
%                   GWAS hits, drug targets) are appended to the gene set.
%                   Default: [] (skipped).
%   penalty_scale - Safety multiplier applied to the penalty base:
%                     auto_penalty=true  → P = penalty_scale * norm(Q1,2)
%                                          Default: 2.
%                     auto_penalty=false → P = penalty_scale * max(|Q1|)
%                                          Default: 10 (legacy).
%   auto_penalty  - Logical flag (default: true).
%                   true  = spectral-norm-derived P (recommended).
%                   false = legacy elementwise-max P.
%
% OUTPUTS:
%   Q             - Final QUBO matrix with penalty, G×G (dense, symmetric).
%   Q1            - Base matrix before penalty, G×G.
%   P             - Penalty coefficient used.
%
% AUTHOR: Selim Romero, Texas A&M University

    % --- Default arguments -------------------------------------------------
    if nargin < 6 || isempty(Xnet_target)
        Xnet_target = [];
    end
    if nargin < 7 || isempty(penalty_scale)
        penalty_scale = 2;           % default for spectral mode
    end
    if nargin < 8 || isempty(auto_penalty)
        auto_penalty = true;
    end

    ng = size(Xdiff, 1);

    % --- Validate dimensions -----------------------------------------------
    assert(numel(Vdiff)       == ng, 'Vdiff length must equal number of genes.');
    assert(size(Xmnn_A, 1)   == ng, 'Xmnn_A row dimension mismatch.');
    assert(size(Xmnn_B, 1)   == ng, 'Xmnn_B row dimension mismatch.');
    if ~isempty(Xnet_target)
        assert(size(Xnet_target,1) == ng, 'Xnet_target dimension mismatch.');
    end

    % --- Base QUBO matrix Q1 -----------------------------------------------
    Q1 = Xdiff + diag(Vdiff(:)) - (full(Xmnn_A) + full(Xmnn_B));

    if ~isempty(Xnet_target)
        Q1 = Q1 + Xnet_target;
        fprintf('  Pathway prior X_net: included\n');
    else
        fprintf('  Pathway prior X_net: skipped (pure pathway-gene analysis)\n');
    end

    % --- Cardinality penalty -----------------------------------------------
    if auto_penalty
        % Spectral norm = largest singular value = largest |eigenvalue|
        % for symmetric Q1.  Provides a principled lower bound on P.
        sigma_max = norm(Q1, 2);
        P = penalty_scale * sigma_max;
        fprintf('  Penalty (spectral): ||Q1||_2 = %.4f  →  P = %.4f  (scale = %.1f)\n', ...
                sigma_max, P, penalty_scale);
    else
        P = penalty_scale * max(abs(Q1(:)));
        fprintf('  Penalty (legacy max): P = %.4f  (scale = %.1f)\n', P, penalty_scale);
    end

    % Build penalty matrix efficiently (no for-loops):
    %   all elements = P, then fix diagonal
    Q_penalty = P * ones(ng, ng);
    diag_fix  = P * (1 - 2*K) - P;           % diagonal should be P*(1-2K), not P
    Q_penalty = Q_penalty + diag(diag_fix * ones(ng, 1));

    % --- Final QUBO matrix -------------------------------------------------
    Q = Q1 + Q_penalty;

    % Symmetrise (numerical safety)
    Q = (Q + Q') / 2;

    fprintf('QUBO assembled: %d×%d | K=%d | P=%.4f\n', ng, ng, K, P);
    fprintf('  Q1 range: [%.4f, %.4f]\n', min(Q1(:)), max(Q1(:)));

end
