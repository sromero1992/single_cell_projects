function [CM, mass_zero] = center_mass(X)
    %{
    center_mass computes the centroid of the data as mj*sum_j( R - x_ij)
    INPUT:
    X ------> Count matrix in a gene by cell fashion
    OUTPUT:
    CM --------> Center of mass for input count matrix
    mass_zero -> Boolean for loop control
    %}
    X = normr(X);
    num_cells = size(X,1);
    num_genes = size(X,2);
    % Mass must be the sum of all genes in each cell (cell mass)
    % Sequencing depth for each cell
    MM = sum(X, 2);
    % All mass or all the sequencing depths
    M = sum(MM);

    % Handle zero-mass case upfront
    mass_zero = M == 0;
    if mass_zero
        CM = zeros(num_genes, num_cells);  % No calculations needed if mass is zero
        return;
    end

    % Getting center of mass for a subsample set of cells? (here are all cells)
    % Center of mass (vector) RI = sum_j^ncells mj * x_ij) for i-th gene 
    R = MM'*X/M; % <R> expectation of "position"-gene expression

    % Center of mass equation must be sum_j( RI - mj * x_ij)
    % Calculate CM matrix efficiently using broadcasting
    CM = bsxfun(@minus, R, X) .* MM;

end