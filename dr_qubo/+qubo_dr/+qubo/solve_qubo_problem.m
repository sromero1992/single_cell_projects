function [selected_idx, sub_g_net, result] = solve_qubo_problem(Q, g)
% SOLVE_QUBO_PROBLEM Solve the QUBO problem using MATLAB's optimization toolbox.
%
%   [selected_idx, sub_g_net, result] = SOLVE_QUBO_PROBLEM(Q, g)
%
% Solves the quadratic unconstrained binary optimization (QUBO) problem
% using MATLAB's qubo() and solve() functions.
%
% INPUT:
%   Q   - QUBO matrix (ng × ng, symmetric)
%   g   - Gene names (string array, length ng)
%
% OUTPUT:
%   selected_idx  - Boolean index of selected genes (length ng)
%   sub_g_net     - Names of selected genes (string array)
%   result        - Optimization result object
%
% AUTHOR: Selim Romero, Texas A&M University

    g = string(g(:));
    ng = length(g);

    % Create QUBO problem
    qprob = qubo(Q);

    % Solve the problem
    % This uses MATLAB's built-in QUBO solver (Optimization Toolbox)
    options = optimoptions('solve', 'Display', 'iter', 'MaxTime', 300);
    result = solve(qprob);

    % Extract solution
    x_sol = result.BestX;
    selected_idx = (x_sol == 1);

    % Get selected gene names
    sub_g_net = g(selected_idx);

    % Print summary
    num_selected = nnz(selected_idx);
    fprintf('QUBO solution found:\n');
    fprintf('  Objective value: %.6f\n', result.BestFval);
    fprintf('  Number of selected genes: %d\n', num_selected);
    fprintf('  Selected genes: %s\n', strjoin(sub_g_net, ', '));

end
