% Extract Expression Data
%
% This script processes single-cell RNA-sequencing data to analyze the gene 
% expression of a specified cell type (e.g., Macrophage). It loads data 
% from a .mat file and an Excel file, filters the relevant expression 
% data, normalizes it, performs clustering if specified, and generates a 
% heatmap of mean gene expression.
%
% Inputs:
% - "Bulk_Analysis.mat": A .mat file containing a SingleCellExperiment object.
% - An Excel file containing gene expression data across different sheets for various cell types.
%
% Outputs:
% - A heatmap visualizing the mean gene expression for the specified cell type 
%   across different time batches.
% - A CSV file with clustering results.

% Load bulk analysis data
% data = load("Bulk_Analysis.mat");   % Load SingleCellExperiment data
% sce = data.sce;                     % Extract SingleCellExperiment object
% clear data;

% Specify the cell type and load its specific analysis file
celltype = "T cells";  % Specify the cell type of interest
fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
D = readtable(fname, 'ReadVariableNames', true);

% Sort rows based on acrophase and amplitude
D = sortrows(D, ["Acrophase_24" "Abs_Amp"], ["ascend" "descend"]);

% Filter adjusted p-values and acrophase threshold
Dwork = D(D.pvalue <= 0.05, :);      % Filter for significant adjusted p-values
%Dwork = D(D.pvalue_adj <= 0.05, :);      % Filter for significant adjusted p-values
Dwork = Dwork(Dwork.Amp >= 0, :);        % Filter for acrophase threshold
%Dwork = Dwork(Dwork.Amp < 0, :);        % Filter for acrophase threshold

% Define list of classical circadian genes (partial matches allowed)
classic_circ = ["arn" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
    "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
    "nfkb" "per" "wnt" "nrd" "rev" "pik"];

% Filter `Dwork` to include only rows with genes containing any substring in `classic_circ`
glcirc = contains(lower(Dwork.Genes), classic_circ, 'IgnoreCase', true);
%Dwork = Dwork(glcirc, :);

% Extract circadian gene list for the specified cell type
gl = string(Dwork.Genes);                % Convert gene names to a string array
ic = strcmpi(sce.c_cell_type_tx, celltype);  % Identify indices of the specified cell type
ig = contains(lower(sce.g), classic_circ);   % Check which genes in `sce.g` contain circadian genes

% Extract expression data for the circadian genes and specified cell type
X = sce.X(ig, ic);                        % Filter expression matrix for circadian genes
g = sce.g(ig);                            % Get gene names for filtered circadian genes


% Extract gene list for the specified cell type
gl = string(Dwork.Genes);                % Convert gene names to string array
ic = sce.c_cell_type_tx == celltype;     % Identify indices of the specified cell type
ig = ismember(sce.g, gl);                % Check which genes are in the list from the specified sheet

% Extract expression data for the selected genes and cell type
X = sce.X(ig, ic);                       % Filter expression matrix
X = sc_norm(X);
%X = sc_impute(X);
g = sce.g(ig);                           % Get gene names for filtered genes
sce_new = SingleCellExperiment(X, g);    % Create new SingleCellExperiment for filtered data
sce_new.c_batch_id = sce.c_batch_id(ic); % Update batch information
sce_new.c_cell_type_tx = sce.c_cell_type_tx(ic);
sce_new.c_cell_id = sce.c_cell_id(ic);   % Update cell information

% Calculate mean expression across batches
time_batches = unique(sce_new.c_batch_id); % Unique time batches
nt = length(time_batches);                % Number of time batches
ng = length(sce_new.g);                   % Number of genes
gzts = zeros(nt, ng);                     % Initialize mean expression matrix

for igs = 1:ng
    for ib = 1:nt
        ic = sce_new.c_batch_id == time_batches(ib);  % Indices for current batch
        gzts(ib, igs) = mean(sce_new.X(igs, ic),'omitnan');     % Mean expression across batches
    end
end

% Prepare data for analysis
data = gzts';                             % Transpose for (genes x batches)

% Scaling data
norm_type = 'zscore';                      % Choose normalization type
switch norm_type
    case 'zscore'
        data_scaled = normalize(data', "zscore", "std"); % Z-score normalization
    case 'norm' 
        data_scaled = normalize(data', "norm", 2);       % L2 normalization
    case 'scale'
        data_scaled = normalize(data', "range", [-1, 1]); % Min-max scaling
    case 'center' 
        data_scaled = normalize(data', "center", "mean"); % Centering
    otherwise
        data_scaled = normalize(data', "zscore", "std"); % Default to z-score
end
data_scaled = data_scaled';               % Transpose back to (batches x genes)

% Define clustering method
clusteringMethod = 'spectral';  % Clustering method ('spectral', 'hierarchical', 'none')

switch clusteringMethod
    case 'spectral'
        % Spectral clustering parameter search
        num_clusters_range = [5 10 20];       % Range of clusters
        sigma_values = [1 2 3];               % Sigma values for similarity graph
        best_silhouette = -1;                 % Initialize silhouette score
        best_num_clusters = num_clusters_range(1);
        best_sigma = sigma_values(1);

        for num_clusters = num_clusters_range
            for sigma = sigma_values
                similarity_graph = exp(-pdist2(data_scaled, data_scaled).^2 / (2 * sigma^2));
                cluster_idx = spectralcluster(similarity_graph, num_clusters);
                sil_score = mean(silhouette(data_scaled, cluster_idx));
                if sil_score > best_silhouette
                    best_silhouette = sil_score;
                    best_num_clusters = num_clusters;
                    best_sigma = sigma;
                end
            end
        end
        fprintf("Best sigma %f and num clusters %d \n", best_sigma, best_num_clusters);
        similarity_graph = exp(-pdist2(data_scaled, data_scaled).^2 / (2 * best_sigma^2));
        cluster_idx = spectralcluster(similarity_graph, best_num_clusters);
        [~, sort_idx] = sort(cluster_idx);    % Sort rows by clusters

    case 'hierarchical'
        Z = linkage(data_scaled, 'average');  % Hierarchical clustering with average linkage
        cluster_idx = cluster(Z, 'Cutoff', 0.6, 'Criterion', 'distance'); 
        [~, sort_idx] = sort(cluster_idx);    % Sort rows by clusters

    case 'none'
        sort_idx = 1:size(data_scaled, 1);    % No sorting, keep original order
        % Keep acrophase idx
        for ig = 1:height(data_scaled)
            sort_idx(ig) = find(sce_new.g == gl(ig));
        end
        
    otherwise
        error('Invalid clustering method specified.');
end

% Sort data and labels
data_sorted = data_scaled(sort_idx, :);       % Sort data by specified clustering
yLabels = string(sce_new.g);                  % Gene labels before sorting
yLabels_sorted = yLabels(sort_idx);           % Sorted gene labels

% Check alignment of cluster_idx with Dwork
if ismember(clusteringMethod, {'spectral', 'hierarchical'})
    % Update `Dwork.cluster` only if clustering is performed
    if length(sort_idx) == height(Dwork)
        Dwork.cluster(sort_idx) = cluster_idx; % Assign clusters if aligned
    else
        warning('sort_idx length and Dwork rows differ; cluster assignment skipped.');
    end
end

% Create heatmap
figure;
h = imagesc(data_sorted);                     % Plot heatmap

% Set axis labels and customize appearance
axesHandle = get(h, 'Parent');
axesHandle.XTick = 1:length(time_batches);
axesHandle.YTick = 1:length(yLabels_sorted);
axesHandle.XTickLabels = string(time_batches);
axesHandle.YTickLabels = yLabels_sorted;

title('Mean Gene Expression Heatmap for Specified Cell Type'); % Title
xlabel('Time Batches');                                         % X-axis label
ylabel('Genes');                                                % Y-axis label

% Define the custom blue-white-red colormap
blue = [0, 0, 1];       % RGB for blue
white = [1, 1, 1];      % RGB for white
red = [1, 0, 0];        % RGB for red

n = 256; % Number of colors in the colormap for smooth transitions
blueToWhite = [linspace(blue(1), white(1), n/2)', linspace(blue(2), white(2), n/2)', linspace(blue(3), white(3), n/2)'];
whiteToRed = [linspace(white(1), red(1), n/2)', linspace(white(2), red(2), n/2)', linspace(white(3), red(3), n/2)'];

customColormap = [blueToWhite; whiteToRed]; % Combine blue-white and white-red gradients

% Apply the colormap and adjust the color axis
colormap(customColormap);
% Determine the color axis limits based on the range of data_sorted
dataMin = min(data_sorted(:));
dataMax = max(data_sorted(:));
caxis([dataMin, dataMax]); % Set color axis limits to match the data range
colorbar;       % Add a colorbar for reference
set(gca, 'FontSize', 10); % Increase font size for readability
axis tight;     % Adjust axes tightly

% Save output to CSV
output_fname = strcat(celltype, "_with_clusters.csv");
writetable(Dwork, output_fname);

