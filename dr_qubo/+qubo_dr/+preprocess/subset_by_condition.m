function [Xko, Xwt, cs_ko, cs_wt] = subset_by_condition(X, batch_id, cell_state, ko_label, wt_label)
% SUBSET_BY_CONDITION Extract condition-specific data subsets from scRNA-seq matrix.
%
%   [Xko, Xwt, cs_ko, cs_wt] = SUBSET_BY_CONDITION(X, batch_id, cell_state, ko_label, wt_label)
%
% INPUT:
%   X           - Count matrix (genes × cells)
%   batch_id    - Condition label string array per cell (length = num_cells)
%   cell_state  - Cell state vector (length = num_cells), can be zeros
%   ko_label    - String label for KO condition (default "KO")
%   wt_label    - String label for WT condition (default "WT")
%
% OUTPUT:
%   Xko         - KO condition submatrix (genes × ko_cells)
%   Xwt         - WT condition submatrix (genes × wt_cells)
%   cs_ko       - Cell state vector for KO cells
%   cs_wt       - Cell state vector for WT cells
%
% AUTHOR: Selim Romero, Texas A&M University

    % Input validation
    if nargin < 4
        ko_label = "KO";
    end
    if nargin < 5
        wt_label = "WT";
    end

    batch_id = string(batch_id(:));
    cell_state = cell_state(:);

    % Find indices
    idx_ko = contains(batch_id, ko_label);
    idx_wt = contains(batch_id, wt_label);

    % Extract submatrices
    Xko = X(:, idx_ko);
    Xwt = X(:, idx_wt);

    % Extract cell state vectors
    cs_ko = cell_state(idx_ko);
    cs_wt = cell_state(idx_wt);

end
