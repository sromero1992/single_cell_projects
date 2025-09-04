function Jmat = jaccard_all_pairs(X, binarize)
    %{
    jaccard_all_pairs utilizes count matrix from single-cell RNA sequencing
    data to produce cell-pairs matrix. J(A,B) = |A n B| / (A u B)
    INPUT:
    X --------> Count matrix in cell by gene fashion (c,g) 
    binarize -> Boolean to make X binarized with default 0.2 threshold
    OUTPUT: 
    Jmat --> Jaccard matrix containing cell-cell pair (c,c)
    %}
    % Get the number of rows and columns in the data matrix
    if nargin < 2
        fprintf('jaccard_all_pairs non-binary by default');
        binarize = false;
    end
    [nrows, ~] = size(X);

    % Normalization for each cell (row-wise)
    X = normr(X);
    % Binarize the matrix if required
    if binarize
        X = X > 0.1; % Binarize with threshold of 0.1
    end

    % Calculate row and column wise dot products efficiently using bsxfun
    X_dot_row = sum(X.^2, 2); % Efficient row-wise dot product

    % Calculate intersection and union using matrix operations
    % Intersection similar to overlap when non-binary values
    intersection = X * X'; % Efficient matrix multiplication for intersection
    union = bsxfun(@plus, X_dot_row, X_dot_row') - intersection;

    % Calculate Jaccard index and fill the matrix efficiently
    Jmat = intersection ./ union;
end