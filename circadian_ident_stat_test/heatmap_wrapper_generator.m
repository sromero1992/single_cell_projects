% wrapper
% data = load("Bulk_Analysis.mat");
% sce = data.sce;
% clear data;
% 
% celltype = "Macrophage";
% amp_pos = true;
% circ = false;
% clus_method = 'none'; % 'spectral', 'hierarchical', 'none'
% generateHeatmap_circ(sce, celltype, amp_pos, customName, circ, clus_method)

% Define parameters for benchmarking
celltypes = unique(sce.c_cell_type_tx)';
%clus_methods = {'none', 'spectral'};  % Available clustering methods
clus_methods = {'phase_sort'};  % Available clustering methods
strict = [true false];             % Strict (padj) or not pval filter options
circ_bools = [false, true];            % Circadian filter options
customName = "heatmap_process";
csvname = '_period_24__macro_circadian_analysis.csv';

originalDir = pwd;
% Iterate over each cell type, creating a unique folder per cell type
for icelltype = celltypes
    fprintf("Cell type: %s\n", icelltype);

    % Create the cell-type-specific subdirectory once
    celltypeDir = fullfile('results', char(icelltype));
    if ~exist(celltypeDir, 'dir')
        mkdir(celltypeDir);
    end

    % Specify the source gene expression file directly in the "results" folder
    % Define the full path for the CSV file from the main directory
    srcFile = fullfile(originalDir, strcat(icelltype, csvname));

    % Check if the source file exists before proceeding
    if ~isfile(srcFile)
        warning('Source file %s not found.', srcFile);
        continue;  % Skip this cell type if the file is missing
    end

    cd(celltypeDir)

    % Iterate over combinations of parameters within each cell type
    for iclus = clus_methods
        fprintf("  Clustering method: %s\n", iclus{1});

        % Create a subdirectory for each clustering method within the cell type directory
        clusDir = fullfile('./', iclus{1});
        if ~exist(clusDir, 'dir')
            mkdir(clusDir);
        end

        % Copy the gene expression CSV file into the current clustering method directory
        copyfile(srcFile, clusDir);

        % Move into the clustering method directory for saving outputs
        cd(clusDir);

        % Iterate over amplitude and circadian filter combinations
        for istrict = strict
            fprintf("    Strict filter: %d\n", istrict);
            for ibool = circ_bools
                fprintf("    Circadian filter: %d\n", ibool);

                % Generate the heatmap with the specified parameters
                generateHeatmap_circ_simple(sce, icelltype, istrict, ...
                                            customName, ibool);
            end
        end

        % Return to the cell-type directory after processing each clustering method
        cd('..');
    end

    % Return to the main directory after processing the cell type
    cd(originalDir);
end



