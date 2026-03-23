function [Xdiff, Vdiff] = compute_differential(Xwt2norm, Xko2norm, Xwt_norm, Xko_norm, cs_wt, cs_ko)
% COMPUTE_DIFFERENTIAL Compute differential co-expression and cell state vectors.
%
%   [Xdiff, Vdiff] = COMPUTE_DIFFERENTIAL(Xwt2norm, Xko2norm, Xwt_norm, Xko_norm, cs_wt, cs_ko)
%
% Computes the differential gene similarity matrix (WT - KO) and the
% differential cell-state-weighted gene expression vector.
%
% INPUT:
%   Xwt2norm    - WT Gram matrix (genes × genes)
%   Xko2norm    - KO Gram matrix (genes × genes)
%   Xwt_norm    - WT normalized expression (genes × wt_cells)
%   Xko_norm    - KO normalized expression (genes × ko_cells)
%   cs_wt       - WT cell state vector (length = wt_cells), can be zeros
%   cs_ko       - KO cell state vector (length = ko_cells), can be zeros
%
% OUTPUT:
%   Xdiff       - Differential matrix (genes × genes), WT - KO
%   Vdiff       - Differential cell-state vector (genes × 1)
%
% NOTE: Negative values in Xdiff indicate increased co-expression in KO.
%
% AUTHOR: Selim Romero, Texas A&M University

    % Compute differential matrix
    Xdiff = Xwt2norm - Xko2norm;

    % Normalize cell state vectors
    cs_wt = cs_wt(:);
    cs_ko = cs_ko(:);

    cs_wt_norm = cs_wt / (norm(cs_wt) + 1e-10);
    cs_ko_norm = cs_ko / (norm(cs_ko) + 1e-10);

    % Compute differential cell-state vector
    % Project normalized expression onto normalized cell state
    v_wt = Xwt_norm * cs_wt_norm;
    v_ko = Xko_norm * cs_ko_norm;

    Vdiff = v_wt - v_ko;

end
