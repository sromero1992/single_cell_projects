function [p_adj] = bh_adjust_pvalues(pvals)
    % Function to compute Benjamini-Hochberg (BH) adjusted p-values
    % INPUT:
    % pvals => Vector of p-values from F-test
    % OUTPUT:
    % p_adj => Vector of adjusted p-values

    % Number of tests
    m = length(pvals);

    % Sort p-values in ascending order and keep track of original indices
    [sorted_pvals, sort_idx] = sort(pvals);

    % Initialize adjusted p-values
    p_adj_sorted = zeros(m, 1);

    % Compute adjusted p-values using BH procedure
    for i = 1:m
        % Compute the BH adjusted p-value (divide by rank order)
        p_adj_sorted(i) = sorted_pvals(i) * m / i;
    end

    % Ensure monotonicity (p_adj_sorted must be non-decreasing)
    p_adj_sorted = min(cummin(p_adj_sorted(end:-1:1)), 1);
    p_adj_sorted = p_adj_sorted(end:-1:1);

    % Map adjusted p-values back to original pvals order
    p_adj = zeros(m, 1);
    p_adj(sort_idx) = p_adj_sorted;
end
