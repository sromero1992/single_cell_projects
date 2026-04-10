function [dS, Vdiff] = compute_differential(Xref2norm, Xtest2norm, Xref_norm, Xtest_norm, cs_ref, cs_test)
% COMPUTE_DIFFERENTIAL Compute differential co-expression and cell state vectors.
%
%   [dS, Vdiff] = COMPUTE_DIFFERENTIAL(Xref2norm, Xtest2norm, Xref_norm, Xtest_norm, cs_ref, cs_test)
%
% Computes the differential gene similarity matrix (reference − test) and the
% differential cell-state-weighted gene expression vector.
%
% INPUT:
%   Xref2norm   - Reference Gram matrix (genes × genes)
%   Xtest2norm  - Test Gram matrix (genes × genes)
%   Xref_norm   - Reference normalized expression (genes × ref_cells)
%   Xtest_norm  - Test normalized expression (genes × test_cells)
%   cs_ref      - Reference cell state vector (length = ref_cells), can be zeros
%   cs_test     - Test cell state vector (length = test_cells), can be zeros
%
% OUTPUT:
%   dS       - Differential matrix (genes × genes), reference − test.
%                 Positive values: reference co-expression gain.
%                 Negative values: test/disease co-expression gain.
%   Vdiff       - Differential cell-state vector (genes × 1)
%
% NOTE: When plotting, negate dS (i.e. plot −dS) so that positive
%       values (red) represent test/disease co-expression gain — the
%       biologically natural direction.
%
% AUTHOR: Selim Romero, Texas A&M University

    % Compute differential matrix: reference − test
    dS = Xref2norm - Xtest2norm;

    % Normalize cell state vectors
    cs_ref  = cs_ref(:);
    cs_test = cs_test(:);

    cs_ref_norm  = cs_ref  / (norm(cs_ref)  + 1e-10);
    cs_test_norm = cs_test / (norm(cs_test) + 1e-10);

    % Compute differential cell-state vector
    % Project normalized expression onto normalized cell state
    v_ref  = Xref_norm  * cs_ref_norm;
    v_test = Xtest_norm * cs_test_norm;

    Vdiff = v_ref - v_test;

end
