 % Include the heatmap generation function here (or ensure it's in your MATLAB path)
    function generateHeatmap_circ_simple(celltype, strict, customName, circ)

        % Specify the file name for gene expression data
        fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
        fname2 = strcat(celltype, "_period_24__macro_circadian_ZTs.csv");

        % Check if the analysis file exists before trying to read it
        if ~exist(fname, 'file') | ~exist(fname2,'file');
            warning('Analysis file "%s" not found. Cannot generate heatmap.', fname);
            return; % Exit the function if the file doesn't exist
        end
        
        D = readtable(fname, 'ReadVariableNames', true);
        Dzts = readtable(fname2, 'ReadVariableNames', true);

        % Sort rows based on acrophase and amplitude
        [D, idx]= sortrows(D, ["Acrophase_24", "Abs_Amp"], ["ascend", "descend"]);
        Dzts = Dzts(idx, :);

        % Filter adjusted p-values and amplitude
        if strict
            % Only F-test
            %idx = and(D.pvalue_corr < 1.0, D.pvalue < 0.05);
            % F-test + Corr t-test
            idx = and(D.pvalue_corr < 0.05, D.pvalue < 0.05);
            Dwork = D(idx, :);
            Dzts = Dzts(idx, :);
        else
            Dwork = D;
        end

        idx = Dwork.Amp < 0;
        Dwork{idx, "Acrophase" } = Dwork{idx, "Acrophase" } + 12;
        Dwork{idx, "Acrophase_24" } = Dwork{idx, "Acrophase_24" } + 12;
        Dwork{idx, "Amp" } = -1*Dwork{idx, "Amp" } + 12;

        Dwork = Dwork(Dwork.Amp >= 0, :);
        Dzts = Dzts(Dwork.Amp >= 0, :);

        % Filter for circadian genes if specified
        if circ
            classic_circ = ["arnt" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
                    "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
                    "nfkb" "per" "wnt" "nr1d" "rev" "pik" "ror"];
            glcirc = contains(lower(Dwork.Genes), classic_circ, 'IgnoreCase', true);
            Dwork = Dwork(glcirc, :);
            Dzts = Dzts(glcirc, :);
        end
        if height(Dwork) <= 1
            disp('Small number of predicted circadian genes...')
            disp(height(Dwork))
            return;
        end

        % % Normalize full set
        % %X = full(sce.X);
        % if strcmp(norm_str, 'lib_size')
        %     % This is regular cells pipeline
        %     %X = sc_norm(full(sce.X));
        %     X = pkg.norm_libsize(sce.X, 1e4);
        %     X = log1p(sce.X);
        % else % 'magic_impute'
        %     % This for cancer cells
        %     X = sc_impute(sce.X, 'MAGIC');
        % end
        % X = sparse(X);
        % sce.X = X;
        % clear X;
        % 
        % % Extract circadian gene list for specified cell type
        % gl = string(Dwork.Genes);
        % ic = strcmpi(sce.c_cell_type_tx, celltype);
        % ig = ismember(lower(sce.g), lower(gl));
        % % Extract expression data for circadian genes and specified cell type
        % X = sce.X(ig, ic);
        % g = sce.g(ig);
        % 
        % sce_new = SingleCellExperiment(X, g);
        % sce_new.c_batch_id = sce.c_batch_id(ic);
        % sce_new.c_cell_type_tx = sce.c_cell_type_tx(ic);
        % sce_new.c_cell_id = sce.c_cell_id(ic);
        % % Calculate mean expression across batches
        % time_batches = unique(sce_new.c_batch_id);
        % nt = length(time_batches);
        % ng = length(sce_new.g);
        % gzts = zeros(nt, ng);
        % for igs = 1:ng
        %     for ib = 1:nt
        %         % Ensure batch ID comparison is robust to type (char/string)
        %         ic = ismember(sce_new.c_batch_id, time_batches(ib));
        %         gzts(ib, igs) = mean(sce_new.X(igs, ic), 'omitnan');
        %     end
        % end
        
        time_batches = string( Dzts.Properties.VariableNames);
        time_batches = time_batches(2:9);
        gzts = Dzts{:, 2:9};
        % Scale data with chosen normalization method
        norm_type = 'zscore';  % Choose normalization type: 'zscore', 'norm', 'scale', 'center'
        switch norm_type
            case 'zscore'
                data_scaled = normalize(gzts', "zscore", "std");
            case 'norm'
                data_scaled = normalize(gzts', "norm", 2);
            case 'scale'
                data_scaled = normalize(gzts', "range", [-1, 1]);
            case 'center'
                data_scaled = normalize(gzts', "center", "mean");
            otherwise
                data_scaled = normalize(gzts', "zscore", "std");
        end
        data_scaled = data_scaled';

        % Get index order based on final sort by "Acrophase_24"
        % % We need to match the genes in Dwork with the genes in sce_new
        % [~, sort_idx] = ismember(lower(Dwork.Genes), lower(sce_new.g));
        % 
        % % Remove indices that were not found (where ismember returned 0)
        % found_genes_idx = sort_idx > 0;
        % sort_idx = sort_idx(found_genes_idx);
        % sorted_gene_names = Dwork.Genes(found_genes_idx); % Keep track of the names
        % 
        % if isempty(sort_idx)
        %     disp('No matching genes found in SCE object after filtering.');
        %     return;
        % end
        % 
        % % Use sort_idx to extract data_sorted and corresponding gene names
        % data_sorted = data_scaled(sort_idx, :);
        sorted_gene_names_sce = string(Dzts{:,1}); % sce_new.g(sort_idx); % Get the gene names from sce_new based on sort_idx

        % Generate heatmap
        figureHandle = figure;
        h = imagesc(data_scaled);
        set(h, 'AlphaData', ~isnan(data_scaled));  % Handle NaN values gracefully

        % Set axis labels and customize appearance
        axesHandle = get(h, 'Parent');
        axesHandle.XTick = 1:length(time_batches);
        axesHandle.XTickLabels = time_batches;
        axesHandle.YTick = 1:length(sorted_gene_names_sce); % Set YTicks for each gene
        axesHandle.YTickLabels = sorted_gene_names_sce;    % Set YTickLabels to sorted gene names
        
        % Set axis labels and customize appearance
        if size(data_scaled, 1) > 20
            axesHandle.YTick = [];
            axesHandle.YTickLabels = [];
        end
        
        title(sprintf('Mean Gene Expression Heatmap for %s', celltype));
        xlabel('Time Batches');
        ylabel('Genes (Sorted by Acrophase)');

        % Define custom blue-white-red colormap
        n = 256;
        blueToWhite = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', linspace(1, 1, n/2)'];
        whiteToRed = [linspace(1, 1, n/2)', linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
        customColormap = [blueToWhite; whiteToRed];

        colormap(customColormap);
        % Set caxis based on the actual data range
        data_range = [min(data_scaled(:)), max(data_scaled(:))];
        if diff(data_range) > 0 % Avoid error if all values are the same
             caxis(data_range);
        end
        colorbar;
        set(gca, 'FontSize', 10);
        axis tight;
        % Construct the filename based on `customName`, `clus_method`, `strict`, and `circ`
        fname_prefix = strcat(celltype);
        if ~isempty(customName)
            fname_prefix = strcat(fname_prefix, '_', customName);
        end
        fname_prefix = strcat(fname_prefix, '_');

        if strict
            fname_prefix = strcat(fname_prefix, '_strict');
        else
            fname_prefix = strcat(fname_prefix, '_not_strict');
        end

        if circ
            fname_prefix = strcat(fname_prefix, '_core_circ');
        end

        % Specify the output file name for the table (using Dwork which is already filtered)
        output_table_fname = strcat(fname_prefix, "_Dwork_analysis_results.csv");
        writetable(Dwork, output_table_fname);

        % Save figure
        saveas(figureHandle, strcat(fname_prefix, '.svg'));
        close(figureHandle);
    end