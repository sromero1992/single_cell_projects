function [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, rm_low_conf, period12, ...
                                    custom_genelist, custom_celltype, plothist)

    % USAGE:
    % times = [0 3 6 9 12 15 18 21]'; 
    % old_labels = unique(sce.c_batch_id);
    % % Convert times to ZT labels
    % new_labels = cell(numel(times), 1);
    % for i = 1:numel(times)
    %     if times(i) >= 0
    %         new_labels{i} = sprintf('ZT%02d', times(i));
    %     else
    %         new_labels{i} = sprintf('%d', times(i));
    %     end
    % end
    % new_labels = string(new_labels);
    % tmeta = table(old_labels, new_labels, times);
    % [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta)
    tic;
    rng('default');
    
    if nargin < 3 || isempty(rm_low_conf); rm_low_conf = true; end
    if nargin < 4 || isempty(period12); period12 = false; end
    if nargin < 5 || isempty(custom_genelist); custom_genelist = {}; end
    if nargin < 6 || isempty(custom_celltype); custom_celltype = {}; end
    if nargin < 7 || isempty(plothist); plothist = true; end

    if period12 
        disp("Circadian identification with 12 hrs period...")
    else
        disp("Circadian identification with 24 hrs period...")
    end

    % Merge replicates for this analysis or just re-label
    batches = unique(sce.c_batch_id);

    tmeta.times = sortrows(tmeta.times);
    %time_cycle = max(tmeta.times);

    % Needs to be generalized... % Based on Sato data
    time_step = mean(diff(tmeta.times)); 
    % Rename batches 
    for ib = 1:length(batches)
        str_idx = find(batches(ib) == tmeta.old_labels);
        idx = find(sce.c_batch_id == batches(ib));
        sce.c_batch_id(idx) = tmeta.new_labels(str_idx);
    end
    batches = unique(sce.c_batch_id);
    
    disp("New batches")
    disp(batches')

    % % Initialize parallel pool if not already running
    % if isempty(gcp('nocreate'))
    %     numCores = ceil(feature('numcores')/2);
    %     numCores = max(2,numCores);
    %     parpool(numCores);
    %     disp(['Parallel pool initialized with ', num2str(numCores), ' cores.']);
    % 
    %     % Warm up the parallel pool (optional but recommended)
    %     disp('Warming up parallel pool...');
    %     tic; % Start timer for warm-up
    %     parfor i = 1:numCores
    %         pause(0.01); % Small task to initialize workers
    %     end
    %     warmUpTime = toc; % Stop timer and get elapsed time
    %     disp(['Parallel pool warmed up in ', num2str(warmUpTime), ' seconds.']);
    % end
    
    % All cell types available
    cell_type_list = unique(sce.c_cell_type_tx);
    ncell_types = length(cell_type_list);

    % Number of time points
    nztps = length(unique(tmeta.new_labels));
    
    info_p_type = zeros(ncell_types, 6);
    for icell_type = 1:ncell_types   
        % Extract count matrix for ith cell type
        cell_type = cell_type_list(icell_type);
        if ~isempty(custom_celltype)
            if ~ismember(cell_type, custom_celltype); continue; end
        end

        fprintf("Processing cell type %s \n", cell_type);

        % Gene genes and cells information
        idx = find(sce.c_cell_type_tx == cell_type);
        sce_sub = sce.selectcells(idx);
        %sce_sub = sce_sub.qcfilter;
        % Light QC to remove non needed genes and poor cells
        sce_sub = sce_sub.qcfilter(500, 0.20, 10);

        % Normalizing count data for that cell type
        X = full(sce_sub.X);
        X = sc_norm(X);
        X = log1p(X);
        %X = sc_impute(X);
        X = sparse(X);
        sce_sub.X = X;
        clear X idx;

        if isempty(custom_genelist)
            disp("Circadian analysis for all genes");
            gene_list = sce_sub.g;
        else
            disp("Circadian analysis for custom genes");
            % Ensure custom_genelist is a cell array of strings
            if ischar(custom_genelist) || isstring(custom_genelist)
                custom_genelist = cellstr(custom_genelist);
            end
            % Ensure custom_genelist is a row vector
            if iscolumn(custom_genelist)
                custom_genelist = custom_genelist';
            end
            % Match the custom genes with those in sce_sub.g
            [lic, ~] = ismember(custom_genelist, sce_sub.g);
            % Filter only matching genes
            gene_list = custom_genelist(lic);
            % Display matching genes
            disp("Matching genes: " + strjoin(gene_list, ', '));
        end
        ngene = length(gene_list);

        % Batch name for each time point 
        batch_time = unique(sce_sub.c_batch_id);
        % Number of time points
        nzts = length(batch_time);
        fprintf("Number of NZTS after sub-sample: %d \n", nzts);
        
        % If number of time points do not match experimental time points, dump cell type
        if nzts ~= nztps
            disp("Number of time points does not match to input metadata table")
            continue; 
        end
   
        % Gene blocking parameters
        if sce.NumCells > 15000
            block_size = 50; % Adjust block size as needed
        else
            block_size = 500; % Adjust block size as needed
        end
        num_genes = length(gene_list);
        num_blocks = ceil(num_genes / block_size);
    
        progressbar('Starting circadian analysis...'); % Initialize progress bar
        
        %Initialize final arrays.
        acro = zeros(num_genes,1);
        amp = zeros(num_genes,1);
        T = zeros(num_genes,1);
        mesor = zeros(num_genes,1);
        R0 = zeros(num_genes,nzts);
        p_value = zeros(num_genes,1);

        for block_idx = 1:num_blocks
            start_idx = (block_idx - 1) * block_size + 1;
            end_idx = min(block_idx * block_size, num_genes);
            current_gene_block = gene_list(start_idx:end_idx);
            num_genes_in_block = length(current_gene_block);
    
            % Initialize temporary arrays for the current block
            tmp_acro = zeros(num_genes_in_block, 1);
            tmp_amp = zeros(num_genes_in_block, 1);
            tmp_T = zeros(num_genes_in_block, 1);
            tmp_p_value = zeros(num_genes_in_block, 1);
            tmp_R0 = zeros(num_genes_in_block, nzts);
            tmp_mesor = zeros(num_genes_in_block, 1);
    
            % Create gene index lookup table for current block
            gene_indices = zeros(num_genes_in_block,1);
            for gene_index_block = 1:num_genes_in_block
                gene_indices(gene_index_block) = find(sce_sub.g == string(current_gene_block{gene_index_block}));
            end
    
            %parfor igene_block = 1:num_genes_in_block
            for igene_block = 1:num_genes_in_block
                ig = gene_indices(igene_block);
                Xg_zts = {};
                for it = 1:nzts
                    ics = find(sce_sub.c_batch_id == batch_time(it));
                    if ~isempty(ig)
                        Xg_zts{it} = full(sce_sub.X(ig, ics));
                        tmp_R0(igene_block, it) = mean(Xg_zts{it});
                    else
                        Xg_zts{it} = zeros(size(ics));
                        tmp_R0(igene_block, it) = NaN;
                    end
                end
                [tmp_acro(igene_block), tmp_amp(igene_block), tmp_T(igene_block), tmp_mesor(igene_block), tmp_p_value(igene_block)] = ...
                           estimate_phaseR(Xg_zts, time_step, period12, 'Ftest');
            end
    
            % Aggregate results for the current block
            acro(start_idx:end_idx) = tmp_acro;
            amp(start_idx:end_idx) = tmp_amp;
            T(start_idx:end_idx) = tmp_T;
            mesor(start_idx:end_idx) = tmp_mesor;
            R0(start_idx:end_idx, :) = tmp_R0;
            p_value(start_idx:end_idx) = tmp_p_value;
    
            % Clear temporary variables to free memory
            clear tmp_acro tmp_amp tmp_T tmp_mesor tmp_R0 tmp_p_value current_gene_block gene_indices;
    
            % Update progress bar
            progress = (block_idx / num_blocks) * 100;
            progressbar(progress);
        end
        

        clear tmp_mesor tmp_T tmp_amp tmp_acro;

        acro_formatted = acro;
        acro_formatted(acro < 0) = acro(acro < 0) + 24;
        acro_formatted(acro > 24) = acro(acro > 24) - 24;
    
        % Pvalue adjusted with Benjamini-Hochberg (BH) Procedure
        p_adj = bh_adjust_pvalues(p_value);

        T1 = table(gene_list, amp, abs(amp), mesor, acro, ...
                                      acro_formatted, T, p_value, p_adj);
        
        T1.Properties.VariableNames = ["Genes", "Amp", "Abs_Amp", "Mesor", "Acrophase", ...
                                       "Acrophase_24", "Period", "pvalue","pvalue_adj"];
        
        T2 = table(gene_list, R0(:,1), R0(:,2), R0(:,3), R0(:,4), ...
                               R0(:,5), R0(:,6), R0(:,7), R0(:,8));
        
        T2.Properties.VariableNames = ["Genes","ZT00","ZT03","ZT06","ZT09","ZT12","ZT15","ZT18","ZT21"];  

        % Remove NaNs (not optional)
        rm_nans_idx = ~isnan(T1.pvalue);
        T1 = T1(rm_nans_idx, :);
        T2 = T2(rm_nans_idx, :);

        % Remove low confidence genes?
        rm_nconfs_idx = T1.pvalue < 0.05;
        num_adj_conf_g = sum( T1.pvalue_adj <0.05);
        num_conf_g = sum(rm_nconfs_idx);
        num_n_conf_g = sum(~rm_nconfs_idx);
        if rm_low_conf
            T1 = T1(rm_nconfs_idx, :);
            T2 = T2(rm_nconfs_idx, :);
        end
        
        % Sort by acrophase and amplitude to classify better
        [T1, idx] = sortrows(T1, ["pvalue_adj", "pvalue", "Acrophase_24", "Abs_Amp"], ...
                                {'ascend', 'ascend', 'ascend', 'descend'});
        
        T2 = T2(idx, :);
        if period12
            per_label = "_period_12_";
        else
            per_label = "_period_24_";
        end
    
        fname = strcat(cell_type, per_label);
        ftable_name = strcat(fname, "_macro_circadian_analysis.csv");
        writetable(T1, ftable_name);
        
        ftable_name = strcat(fname, "_macro_circadian_ZTs.csv");
        writetable(T2, ftable_name);
    
        if plothist
            generateHeatmap_circ_simple(sce_sub, cell_type, true, "", false)
        end

        info_p_type(icell_type, :) = [sce_sub.NumCells, sce_sub.NumGenes, ...
                                       length(T1.Genes), num_conf_g, ...
                                       num_n_conf_g, num_adj_conf_g];
        
    end

    T0 = table(cell_type_list, info_p_type(:,1), info_p_type(:,2), ...
                info_p_type(:,3), info_p_type(:,4), info_p_type(:,5), info_p_type(:, 6));
    T0.Properties.VariableNames = ["CellType", "NumCells", "NumGenes", "NumCircadian", ...
                                    "NumConfident", "NumNonConfident", "NumAdjConfident"];
    if isempty(custom_celltype) 
        if ncell_types == 1
            custom_celltype = cell_type_list; 
        else
            custom_celltype = 'non_specific';
        end
    end
    fname = strcat(custom_celltype, "_summary_results.csv");
    writetable(T0, fname);

    disp("Processing complete. Total time elapsed: " + num2str(toc) + " seconds");
end