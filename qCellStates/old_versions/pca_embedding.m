X = sce.X;
X = transpose( full(X) );
X_old = X;
X = pearson_residuals(X);
%X = jaccard_all_pairs(X,true);
%X = center
ndim = 2;
my_seed = 1024;
embedding_methods = 'tsne';

if embedding_methods == 'umap'
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


figure('Position', [100 100 800 600]);
hold on;

% You can use scatter plots to visualize embedding techniques 
subplot(2,2,1)
if ndim == 2
    gscatter(s(:, 1), s(:, 2),sce.c_cell_type_tx);
else
    scatter3(s(:, 1), s(:, 2), s(:,3));
end

subplot(2,2,2)
colormap('hot')
imagesc(X)
colorbar
