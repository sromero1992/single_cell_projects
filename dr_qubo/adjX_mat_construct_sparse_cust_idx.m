function adjX = adjX_mat_construct_sparse_cust_idx(X_orig, method, K, cust_idx)
    % INPUT:
    % X_orig -----> Original raw count matrix (genes x cells)
    % method -----> Neighbor method ('knn' or 'mnn')
    % K ----------> Number of neighbors
    % cust_idx ---> Logical array (n_genes x 1) indicating the genes of interest.
    % OUTPUT:
    % adjX -------> Sparse adjacency matrix focused on cust_idx and their neighbors.
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
    if nargin < 4 || ~islogical(cust_idx) || length(cust_idx) ~= size(X_orig, 1)
        error('cust_idx must be a logical array matching the number of genes in X_orig.');
    end

    fprintf("Computing adjacency matrix for custom index! \n");
    tic;

    % Assuming X_orig is genes x cells
    % If neighbors for cells, use X_orig'
    % If neighbors for genes, use X_orig
    X = log1p(X_orig); % genes by cells (so we are finding gene-gene neighbors)

    % PCA/SVD for dimensionality reduction
    num_svd_components = 50; % Number of components for SVD
    if size(X, 2) > num_svd_components % Only perform SVD if number of genes is greater than components
        [U, S, V] = svds(X, num_svd_components);
        X_reduced = U * S; % (num_genes x num_svd_comps)
    else
        X_reduced = X; % No reduction if fewer genes than components
    end
    fprintf("Time for preparing data: %f \n", toc);

    n = size(X_reduced, 1); % This is num_genes
    adjX = sparse(n, n); % Initialize the full N x N matrix

    tic;
    % Compute K-NN for all genes (X_reduced)
    % This step must still be global to correctly identify all neighbors
    [Idx_knn, ~] = knnsearch(X_reduced, X_reduced, 'K', K+1, 'Distance', 'euclidean');
    Idx_knn = Idx_knn(:, 2:end); % Remove self-loop

    % Get the indices of the custom genes
    cust_gene_indices = find(cust_idx);
    num_cust_genes = length(cust_gene_indices);

    switch method
        case 'knn'
            % Populate adjacency matrix for KNN
            % We only care about edges originating from or pointing to genes in cust_idx
            for i_idx = 1:num_cust_genes % Iterate only through the selected custom genes
                i = cust_gene_indices(i_idx); % Get the global index of the current custom gene

                % Add outgoing edges from 'i' to its K neighbors
                adjX(i, Idx_knn(i, :)) = 1;

                % Add incoming edges to 'i' from its K neighbors
                % This makes it an undirected graph where A -> B implies B -> A
                % This effectively "symmetrizes" the KNN graph around the custom genes.
                for neighbor_of_i_idx = 1:K
                    j = Idx_knn(i, neighbor_of_i_idx);
                    adjX(j, i) = 1; % Mark j as a neighbor of i
                end
            end
            % Final symmetrization (redundant if loop logic makes it perfectly symmetric, but safe)
            adjX = max(adjX, adjX');

        case 'mnn'
            % Populate adjacency matrix for MNN
            % We only check for MNN relationships involving at least one gene from cust_idx
            % To ensure symmetry and correct MNN for *any* gene that is a neighbor of a custom gene,
            % we need to check both directions.
            for i_idx = 1:num_cust_genes % Iterate only through the selected custom genes
                i = cust_gene_indices(i_idx); % Get the global index of the current custom gene

                for j_idx_in_list = 1:K % Iterate through the K nearest neighbors of 'i'
                    j = Idx_knn(i, j_idx_in_list); % 'j' is one of the K-NNs of 'i'

                    % Check if 'i' is also a K-NN of 'j'
                    if ismember(i, Idx_knn(j, :))
                        adjX(i, j) = 1; % Mark the mutual connection
                        adjX(j, i) = 1; % Ensure symmetry
                    end
                end
            end
            % No extra max(adjX, adjX') needed here as MNN loop ensures symmetry for relevant edges
            % However, if you only wanted connections *between* genes in cust_idx, that's different.
            % The current logic gets all MNN where *at least one* of the pair is a cust_idx gene.
    end

    fprintf("Total time for neighbors: %f \n", toc);
    adjX = sparse(adjX); % Final sparse conversion
end