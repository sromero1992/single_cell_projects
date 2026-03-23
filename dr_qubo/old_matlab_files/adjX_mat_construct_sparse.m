function adjX = adjX_mat_construct_sparse(X_orig, method, K)
    % INPUT:
    % X_orig -----> Original raw count matrix (genes x cells)
    % method -----> Neighbor method ('knn' or 'mnn')
    % K ----------> Number of neighbors
    % use_hvgs ---> Flag to use highly variable genes (needs nhvgs definition)
    % OUTPUT:
    % adjX -------> Sparse adjacency matrix
    % AUTHOR: Selim Romero, Texas A&M University

    % Input validation
    method = lower(method);
    if ~ismember(method, {'knn', 'mnn'})
        error('Method must be either ''knn'' or ''mnn''.');
    end
    if nargin < 3; K = 15; end
    if ~isnumeric(K) || K <= 0 || K ~= round(K)
        error('K must be a positive integer.');
    end

    fprintf("Computing adjacency matrix! \n");
    tic;

    % Neighboars for cells
    %X = log1p(X_orig)'; % cell by genes
    % Neighboars for genes
    X = log1p(X_orig); % genes by cells

    % PCA/SVD for dimensionality reduction
    num_svd_components = 50; % Number of components for SVD
    if size(X, 2) > num_svd_components % Only perform SVD if number of genes is greater than components
        [U, S, V] = svds(X, num_svd_components); % SVD on X (cells x genes)
        X_reduced = U * S; % (cells x num_svd_comps) | genes x num_svd_comps
    else
        X_reduced = X; % No reduction if fewer genes than components
    end
    fprintf("Time for preparing data: %f \n", toc);

    % num cells or num genes
    n = size(X_reduced, 1);
    adjX = sparse(n, n);

    tic;
    % Compute K-NN for all vars using knnsearch
    % knnsearch returns indices of K nearest neighbors
    
    % For 'euclidean' distance
    [Idx_knn, ~] = knnsearch(X_reduced, X_reduced, 'K', K+1, 'Distance', 'euclidean');
    % [Idx_knn, ~] = knnsearch(X_reduced, X_reduced, 'K', K+1, 'Distance', 'cosine');
    % K+1 because knnsearch includes the point itself as its own nearest neighbor (distance 0)
    % Cosine distance is 1 - cosine similarity. Lower distance is higher similarity.

    % Remove self-loop (the first neighbor is always the point itself)
    Idx_knn = Idx_knn(:, 2:end); % Now Idx_knn is n x K

    switch method
        case 'knn'
            % Populate adjacency matrix for KNN
            for i = 1:n
                adjX(i, Idx_knn(i, :)) = 1;
            end
            adjX = max(adjX, adjX');
        case 'mnn'
            % Populate adjacency matrix for MNN
            for i = 1:n
                for j_idx = 1:K
                    j = Idx_knn(i, j_idx); % j is a K-NN of i
                    
                    % Check if i is a K-NN of j
                    % We need to search if 'i' is present in the K-NN list of 'j'
                    if ismember(i, Idx_knn(j, :))
                        adjX(i, j) = 1;
                        adjX(j, i) = 1; % Ensure symmetry
                    end
                end
            end
    end

    fprintf("Total time for neighbors: %f \n", toc);

    % Final adjacency matrix is sparse
    adjX = sparse(adjX);
end