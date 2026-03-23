function adjX = build_adjacency_sparse(X_orig, method, K)
% BUILD_ADJACENCY_SPARSE Construct KNN/MNN adjacency matrix for full gene set.
%
%   adjX = BUILD_ADJACENCY_SPARSE(X_orig, method, K)
%
% Builds a sparse gene-gene adjacency matrix using KNN or mutual nearest neighbors
% (MNN) on the full input matrix. Uses optional SVD reduction for high-dimensional data.
%
% INPUT:
%   X_orig      - Expression matrix (genes × cells)
%   method      - 'knn' or 'mnn' (default 'mnn')
%   K           - Number of neighbors (default 15)
%
% OUTPUT:
%   adjX        - Sparse symmetric adjacency matrix (genes × genes)
%
% AUTHOR: Selim Romero, Texas A&M University

    if nargin < 2
        method = 'mnn';
    end
    if nargin < 3
        K = 15;
    end

    % Apply log1p transformation
    X_log = log1p(X_orig);

    % SVD reduction if needed (reduce to 50 components if more)
    num_components = min(50, min(size(X_log)) - 1);
    if num_components < min(size(X_log))
        [U, S, ~] = svds(X_log, num_components);
        X_reduced = U * S;
    else
        X_reduced = X_log;
    end

    % Compute pairwise distances
    % Using Euclidean distance: ||x - y||^2
    num_genes = size(X_reduced, 1);
    dist_mat = pdist2(X_reduced, X_reduced, 'euclidean');

    % Find K+1 nearest neighbors (including self)
    [~, idx_knn] = sort(dist_mat, 2);
    idx_knn = idx_knn(:, 1:K+1);

    % Build adjacency matrix based on method
    if strcmpi(method, 'knn')
        % KNN: edge from i to j if j is in KNN(i)
        adjX = sparse(num_genes, num_genes);
        for i = 1:num_genes
            neighbors = idx_knn(i, 2:end);  % Exclude self
            adjX(i, neighbors) = 1;
        end
        % Symmetrize
        adjX = adjX | adjX';

    elseif strcmpi(method, 'mnn')
        % MNN: edge from i to j if j is in KNN(i) AND i is in KNN(j)
        adjX = sparse(num_genes, num_genes);
        for i = 1:num_genes
            neighbors_i = idx_knn(i, 2:end);  % Exclude self
            for j = neighbors_i
                neighbors_j = idx_knn(j, 2:end);
                if any(neighbors_j == i)
                    adjX(i, j) = 1;
                end
            end
        end
        % Symmetrize
        adjX = adjX | adjX';
    else
        error('Unknown method: %s. Use ''knn'' or ''mnn''.', method);
    end

    % Ensure symmetric and sparse
    adjX = sparse(adjX + adjX');
    adjX(adjX > 0) = 1;
    adjX = adjX / 2;  % Remove double-counting from symmetrization

end
