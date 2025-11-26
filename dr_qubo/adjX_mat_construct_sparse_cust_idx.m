function adjX = adjX_mat_construct_sparse_cust_idx(Xsub, method, K)
    % INPUT:
    % Xsub -------> Normalized count matrix of subset genes (genes x cells)
    % method -----> Neighbor method ('knn' or 'mnn')
    % K ----------> Number of neighbors
    % OUTPUT:
    % adjX -------> Sparse adjacency matrix for the subset genes in Xsub.
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
    
    fprintf("Computing adjacency matrix for subset genes! \n");
    tic;
    
    % PCA/SVD for dimensionality reduction
    num_svd_components = 50; % Number of components for SVD
    
    % Use Xsub for SVD
    if size(Xsub, 2) > num_svd_components % Only perform SVD if number of cells is greater than components
        [U, S, ~] = svds(Xsub, num_svd_components);
        X_reduced = U * S; % (num_genes x num_svd_comps)
    else
        X_reduced = Xsub; % No reduction if fewer cells than components
    end
    
    fprintf("Time for preparing data: %f \n", toc);
    
    n = size(X_reduced, 1); % This is num_genes in the subset
    adjX = sparse(n, n); % Initialize the full N x N matrix
    tic;
    
    % Compute K-NN for all genes (X_reduced)
    [Idx_knn, ~] = knnsearch(X_reduced, X_reduced, 'K', K+1, 'Distance', 'euclidean');
    Idx_knn = Idx_knn(:, 2:end); % Remove self-loop
    
    switch method
        case 'knn'
            % Populate adjacency matrix for KNN (ensuring symmetry)
            for i = 1:n 
                % Add outgoing edges from 'i' to its K neighbors
                adjX(i, Idx_knn(i, :)) = 1;
                % Add incoming edges to 'i' from its K neighbors
                for neighbor_of_i_idx = 1:K
                    j = Idx_knn(i, neighbor_of_i_idx);
                    adjX(j, i) = 1; 
                end
            end
            % Final symmetrization (safe practice)
            adjX = max(adjX, adjX');
            
        case 'mnn'
            % Populate adjacency matrix for MNN (Mutual Nearest Neighbor)
            % This loop ensures symmetry as both directions are checked and set.
            for i = 1:n 
                for j_idx_in_list = 1:K % Iterate through the K nearest neighbors of 'i'
                    j = Idx_knn(i, j_idx_in_list); % 'j' is one of the K-NNs of 'i'
                    % Check if 'i' is also a K-NN of 'j'
                    if ismember(i, Idx_knn(j, :))
                        adjX(i, j) = 1; % Mark the mutual connection
                        adjX(j, i) = 1; % Ensure symmetry
                    end
                end
            end
    end
    
    fprintf("Total time for neighbors: %f \n", toc);
    adjX = sparse(adjX); % Final sparse conversion
end