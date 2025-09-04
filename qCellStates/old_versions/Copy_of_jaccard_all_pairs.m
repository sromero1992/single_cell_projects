function Jmat = jaccard_all_pairs(X, binarize)
    %{
    jacard_all_pairs utilizes count matrix from single-cell RNA sequencing
    data to produce cell-pairs matrix. J(A,B) = |A n B| / (A u B)
    INPUT:
    X --------> Count matrix in cell by gene fashion
    binarize -> Boolean to make X binarized with default 0.2 threshold
    OUTPUT: 
    Jmat --> Jaccard matrix containing cell-cell pair 
    %}
    % Get the number of rows and columns in the data matrix
    [nrows, ~] = size(X);

    % Binarizing count matrix ??
    bin_thrsh = 0.1; 
    X = normr(X);
    if binarize 
        % This does not keep intensity or magnitude of countings
        for irow = 1:nrows
            jcols = X(irow,:) > bin_thrsh;
            X(irow,jcols) = 1;
            X(irow,~jcols) = 0;
        end
    end 

    % Create an empty matrix to store Jaccard indices
    Jmat = zeros(nrows, nrows);

    % Loop through all pairs of cells using nested loops
    for i = 1:nrows
        xidot = dot( X(i,:), X(i,:) );
        for j = i+1:nrows
            % Calculate the intersection and union of the cells
            intersection = dot( X(i,:), X(j,:));  
            xjdot = dot( X(j,:), X(j,:) );
            union = xidot + xjdot - intersection;
            % Calculate the Jaccard Index
            Jmat(i, j) = intersection / union;

            % Since the Jaccard Index is symmetrical, fill the lower triangular part
            % of the matrix with the same value
            Jmat(j, i) = Jmat(i, j);
        end
    end
end