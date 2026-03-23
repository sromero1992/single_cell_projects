function [sub_Q_net, sub_Qv] = extract_subnetwork(Q0, Q1, selected_idx, result)
% EXTRACT_SUBNETWORK Extract and format subnetwork matrices from QUBO solution.
%
%   [sub_Q_net, sub_Qv] = EXTRACT_SUBNETWORK(Q0, Q1, selected_idx, result)
%
% Extracts the co-expression matrix and contribution vector for the selected
% subnetwork genes. Negates values for visualization (positive = increased KO co-expression).
%
% INPUT:
%   Q0           - Base QUBO matrix or co-expression matrix (ng × ng)
%   Q1           - Alternative QUBO matrix for contribution (ng × ng)
%   selected_idx - Boolean index of selected genes
%   result       - Optimization result with BestX field
%
% OUTPUT:
%   sub_Q_net    - Subnetwork co-expression matrix (negated, ns × ns)
%   sub_Qv       - Gene contribution vector (ns × 1)
%
% AUTHOR: Selim Romero, Texas A&M University

    % Extract subnetwork matrix (negated for visualization)
    % Positive values after negation = increased KO co-expression
    sub_Q_net = -Q0(selected_idx, selected_idx);

    % Compute gene contribution vector
    % sub_Qv = -Q1 * BestX, then extract selected genes
    x_sol = result.BestX;
    contribution = -Q1 * x_sol;
    sub_Qv = contribution(selected_idx);

    % Ensure outputs are properly shaped
    sub_Q_net = full(sub_Q_net);
    sub_Qv = sub_Qv(:);

    fprintf('Subnetwork extracted:\n');
    fprintf('  Subnetwork size: %d × %d\n', size(sub_Q_net, 1), size(sub_Q_net, 2));
    fprintf('  Contribution vector range: [%.6f, %.6f]\n', min(sub_Qv), max(sub_Qv));

end
