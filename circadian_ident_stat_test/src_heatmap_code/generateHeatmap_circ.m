function generateHeatmap_circ(sce, celltype, amp_pos, customName, circ, clus_method, repres_bool)
    % Specify the file name for gene expression data
    fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
    D = readtable(fname, 'ReadVariableNames', true);

    % Sort rows based on acrophase and amplitude
    D = sortrows(D, ["Acrophase_24", "Abs_Amp"], ["ascend", "descend"]);

    % Filter adjusted p-values and amplitude
    if circ 
        Dwork = D(D.pvalue <= 0.05, :);
    else
        Dwork = D(D.pvalue_adj <= 0.05, :);
    end
    if amp_pos
        Dwork = Dwork(Dwork.Amp >= 0, :);
    else
        Dwork = Dwork(Dwork.Amp < 0, :);
    end
    if height( Dwork) <= 1
        disp('Small number of predicted circadian genes...')
        disp( height( Dwork ))
        return;
    end
    % Filter for circadian genes if specified
    if circ 
        classic_circ = ["arn" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
                        "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
                        "nfkb" "per" "wnt" "nrd" "rev" "pik"];
        glcirc = contains(lower(Dwork.Genes), classic_circ, 'IgnoreCase', true);
        Dwork = Dwork(glcirc, :);
    elseif ~ismember(clus_method, {'spectral', 'hierarchical'}) && repres_bool
        % Pre selection of representative genes per phase values
        %Dwork = Dwork(Dwork.Mesor > 0.1, :);
        nsize = height(Dwork);
        nmaxgene = 70;
        fprintf("Reducing summary heatmap from %d to %d \n", nsize, nmaxgene);
        rng(123);
        if nsize > nmaxgene
            % Define time gap and slices
            time_gap = 3;  % each slice covers 3 hours
            time_slice = 24 / time_gap;  % Number of time slices in a 24-hour period
            set_time = floor(nmaxgene / time_slice);  % Number of genes to select per time slice
            
            % Create temporary table to store selected genes
            Dtmp = table();
        
            for itime = 0:time_gap:24
                % Find indices of genes in the current time slice
                idx = find(Dwork.Acrophase_24 > itime & Dwork.Acrophase_24 <= itime + time_gap);
                
                % Check if there are enough genes in this time slice
                if length(idx) > set_time
                    % Randomly permute indices and select the desired number
                    idx = idx(randperm(length(idx), set_time));
                else
                    % If fewer than set_time genes are available, randomly permute and use all of them
                    idx = idx(randperm(length(idx)));
                end
                
                % Append selected genes to Dtmp
                Dtmp = [Dtmp; Dwork(idx, :)];
            end
            
            % Update Dwork to the reduced selection of 50 genes
            Dwork = Dtmp;
        end

    end

    % Extract circadian gene list for specified cell type
    gl = string(Dwork.Genes);
    ic = strcmpi(sce.c_cell_type_tx, celltype);  
    ig = ismember(lower(sce.g), lower(gl));  

    % Extract expression data for circadian genes and specified cell type
    X = sce.X(ig, ic);
    g = sce.g(ig);

    % Normalize data and create new SingleCellExperiment object
    X = sc_norm(X);
    sce_new = SingleCellExperiment(X, g);
    sce_new.c_batch_id = sce.c_batch_id(ic);
    sce_new.c_cell_type_tx = sce.c_cell_type_tx(ic);
    sce_new.c_cell_id = sce.c_cell_id(ic);

    % Calculate mean expression across batches
    time_batches = unique(sce_new.c_batch_id); 
    nt = length(time_batches);
    ng = length(sce_new.g);
    gzts = zeros(nt, ng);

    for igs = 1:ng
        for ib = 1:nt
            ic = sce_new.c_batch_id == time_batches(ib);
            gzts(ib, igs) = mean(sce_new.X(igs, ic), 'omitnan');
        end
    end

    % Scale data with chosen normalization method
    norm_type = 'zscore';   % Choose normalization type: 'zscore', 'norm', 'scale', 'center'
    switch norm_type
        case 'zscore'
            data_scaled = normalize(gzts, "zscore", "std");
        case 'norm'
            data_scaled = normalize(gzts, "norm", 2);
        case 'scale'
            data_scaled = normalize(gzts, "range", [-1, 1]);
        case 'center'
            data_scaled = normalize(gzts, "center", "mean");
        otherwise
            data_scaled = normalize(gzts, "zscore", "std");
    end
    data_scaled = data_scaled';

    % Define a minimum threshold for the number of rows
    min_rows_for_clustering = 30;  % Set your minimum threshold
    
    % Check the number of rows before clustering
    if ismember(clus_method, {'spectral', 'hierarchical'}) && (size(data_scaled, 1) < min_rows_for_clustering )
        clus_method = 'none';  % Fall back to no clustering
        warning('Insufficient data for clustering; using no clustering instead.');
    end

    % Define clustering method
    switch clus_method
        case 'spectral'
            % Spectral clustering parameter search
            num_clusters_range = [5, 10, 20];
            sigma_values = [1, 2, 3];
            best_silhouette = -1;

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
            fprintf("Optimal num clusters %d and sigma %f \n", best_num_clusters, best_sigma);
            similarity_graph = exp(-pdist2(data_scaled, data_scaled).^2 / (2 * best_sigma^2));
            cluster_idx = spectralcluster(similarity_graph, best_num_clusters);
            [~, sort_idx] = sort(cluster_idx);

        case 'hierarchical'
            Z = linkage(data_scaled, 'average');
            cluster_idx = cluster(Z, 'Cutoff', 0.6, 'Criterion', 'distance'); 
            [~, sort_idx] = sort(cluster_idx);

        case 'none'
            sort_idx = 1:size(data_scaled, 1);

        otherwise
            error('Invalid clustering method specified.');
    end

    % Select genes for heatmap display after clustering
    if ismember(clus_method, {'spectral', 'hierarchical'})
        maxgenes = 70;
        Dwork.cluster = cluster_idx;  % Add cluster information
        clust = unique(Dwork.cluster);
        nclus = length(clust);
        gpc = floor(maxgenes / nclus);
    
        Dwork = sortrows(Dwork, ["cluster", "Mesor"], ["ascend", "descend"]);
        selectedGenes = table();
    
        for i = 1:nclus
            clusterData = Dwork(Dwork.cluster == clust(i), :);
            topGenes = clusterData(1:min(gpc, height(clusterData)), :);
            selectedGenes = [selectedGenes; topGenes];
        end
    
        Dwork_old = Dwork;
        % Now set Dwork to the selected genes
        Dwork = selectedGenes;  % Use selected genes for further processing
    else
        Dwork_old = [];
        Dwork = Dwork;  % If no clustering, still use Dwork
    end
    
    % Final sort by "Acrophase_24"
    Dwork = sortrows(Dwork, ["Acrophase_24", "Abs_Amp"], ["ascend", "descend"]);

    % Match sorted indices of selected genes in sce_new.g
    % Sort one last time by acrophase
    sort_idx = NaN(length(Dwork.Genes), 1);  % Initialize sort_idx with NaN
    gl = string(Dwork.Genes);  % Ensure gl contains the gene names from Dwork
    
    for ig = 1:length(gl)
        idx = find(sce_new.g == gl(ig), 1);  % Find the index of the gene
        if ~isempty(idx)  % Check if the gene was found
            sort_idx(ig) = idx;  % Assign index if found
        else
            % Warn if not found
            warning('Gene %s not found in the SingleCellExperiment object.', gl(ig));
            sort_idx(ig) = NaN;  % Assign NaN or handle appropriately
        end
    end
    
    % Remove NaN values from sort_idx
    sort_idx(isnan(sort_idx)) = [];
    
    % Use sort_idx to extract data_sorted and yLabels_sorted
    data_sorted = data_scaled(sort_idx, :);  
    %yLabels_sorted = sce_new.g(sort_idx);
    
    % Generate heatmap
    figureHandle = figure;
    h = imagesc(data_sorted);
    set(h, 'AlphaData', ~isnan(data_sorted));  % Handle NaN values gracefully
    
    % Set axis labels and customize appearance
    axesHandle = get(h, 'Parent');
    axesHandle.XTick = 1:length(time_batches);
    axesHandle.XTickLabels = string(time_batches);
    axesHandle.YTick = [];               % Remove y-axis tick marks
    axesHandle.YTickLabels = [];          % Remove y-axis tick labels
    
    title(sprintf('Mean Gene Expression Heatmap for %s', celltype));
    xlabel('Time Batches');
    ylabel('Gene Number');
    
    % Define custom blue-white-red colormap
    n = 256;
    blueToWhite = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', linspace(1, 1, n/2)'];
    whiteToRed = [linspace(1, 1, n/2)', linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
    customColormap = [blueToWhite; whiteToRed];
    
    colormap(customColormap);
    caxis([min(data_sorted(:)), max(data_sorted(:))]);
    colorbar;
    set(gca, 'FontSize', 10);
    axis tight;

    
    % Construct the filename based on `amp_pos` and `circ` parameters
    customName = strcat(customName, '_', clus_method);
    
    if amp_pos
        fname = strcat(celltype, '_', customName, '_pos_amp_');
    else
        fname = strcat(celltype, '_', customName, '_neg_amp_');
    end
    
    if circ
        fname = strcat(fname, '_core_circ_');
    end
    
    % Specify the output file name for the table
    output_fname = strcat(fname, "_with_clusters.csv");
    writetable(Dwork, output_fname);
    if ~isempty(Dwork_old)
        output_fname_old = strcat(fname, "_with_clusters_old.csv");
        writetable(Dwork_old, output_fname_old);
    end
    
    % Maximize and save figure
    % Get screen size to help with dynamic positioning
    screenSize = get(0, 'ScreenSize');
    
    % Define desired figure dimensions
    % Maximize vertical size but reduce width for a square-like appearance
    figureWidth = screenSize(3) * 0.2;  % 40% of screen width
    figureHeight = screenSize(4) * 0.9; % 90% of screen height for maximized vertical
    
    % Set figure properties
    set(figureHandle, 'Position', [(screenSize(3) - figureWidth) / 2, ... % Center horizontally
                                   screenSize(4) * 0.05, ...              % Small offset from top of screen
                                   figureWidth, figureHeight]);
    
    % Save figure
    saveas(figureHandle, strcat(fname, '.svg'));
    close(figureHandle);
end
