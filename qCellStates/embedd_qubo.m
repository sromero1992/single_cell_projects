function  embedded = embedd_qubo(X, genelist, genes_qubo, form, method, plotit)
    %{
    This ...
    %}

    fprintf("Using %s for %s\n", method, form);

    X = normr(X);
    % Generate matrix X for embedding
    switch method
        case 'cm'
            [X, ~] = center_mass(X);
        case 'pc'
            X = pearson_residuals(X);
        case 'jac'
            X = jaccard_all_pairs(X, false);
        otherwise
            X = pearson_residuals(X);
    end

    switch form
        case 'cells'
            X = full(X)'; % Cell by gene matrix
            %[~, idx, ~] = intersect(genelist, genes_qubo,'stable');
            %X = X(:,idx);
        case 'genes'
            X = full(X);
        otherwise
            fprintf('Non valid form: $s \n', form);
            return;
    end
    fprintf("Size of X counts %d %d \n", size(X));

    ndim = 2;
    my_seed = 1024;
    embedding_methods = 'tsne';
    if matches(embedding_methods,'umap')
        % Perform PCA using scikit-learn library (assuming it's installed)
        %[coeff, score, latent] = pca(X);
        % Select the number of principal components (e.g., choose based on explained variance)
        %n_components = 50; % Choose an appropriate number of components
        % Reduce data to chosen number of components
        %reduced_data = score(:, 1:n_components);
        s = run.mt_UMAP(X',2);
    else
        % 'exact' algorithm optimizes the Kullback-Leibler divergence of
        % distributions between the original space and the embedded space.
        % knn algorithm
        rng(my_seed)
        s = tsne(X, 'NumDimensions', ndim,'Algorithm', 'barneshut', ...
            'NumPCAComponents', 50,'Standardize', false);
    end
    
    switch form
        case 'cells'
            cluster = zeros(size(X,1),1);
        case 'genes'
            % Classified QUBO genes
            [~, idx, ~] = intersect(genelist, genes_qubo,'stable');
            cluster = zeros(size(genelist,1),1);
            cluster(idx) = 1;
    end
            
    fprintf("idx size: %d %d \n",size(cluster));
    fprintf("s size: %d %d \n",size(s));

    if plotit
        figure('Position', [100 100 800 600]);
        hold on;
        
        % You can use scatter plots to visualize embedding techniques
        subplot(2,2,1)
        if ndim == 2
            gscatter(s(:, 1), s(:, 2), cluster);
            %gscatter(s(:, 1), s(:, 2),sce.c_cell_type_tx);
        else
            scatter3(s(:, 1), s(:, 2), s(:,3));
        end
        
        subplot(2,2,2)
        colormap('hot')
        imagesc(X)
        colorbar
    end
end