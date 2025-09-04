function sce_circ_plot_gene(sce, tmeta, cust_cells, period12, cust_gene, ...
                            axHandle, print_scdata, norm_str, use_violin_plot)
    % Plot identified circadian gene for the specified custom gene
    % Inputs:
    %   sce: SingleCellExperiment object
    %   tmeta: Time metadata (struct or table)
    %   cust_cells: Custom cell type to plot
    %   period12: Boolean, if true, use 12-hour period estimation
    %   cust_gene: Custom gene name to plot
    %   axHandle: Axes handle to plot on
    %   print_scdata: Boolean, if true, print single-cell data (scatter or violin)
    %   norm_str: Normalization method ('lib_size' or 'magic_impute')
    %   use_violin_plot: Boolean, if true, plot violin plots; otherwise, plot scatter plots.

    % Ensure required inputs
    if nargin < 4 || isempty(period12)
        period12 = false;
    end
    if nargin < 5 || isempty(cust_gene)
        error('Custom gene must be specified.');
    end
    if nargin < 6 || isempty(axHandle)
        error('Axes handle must be specified.');
    end
    if nargin < 7 || isempty(print_scdata)
        print_scdata = false;
    end
    if nargin < 8 || isempty(norm_str)
        norm_str = 'lib_size';
    end
    if nargin < 9 || isempty(use_violin_plot)
        use_violin_plot = false; % Default to scatter plot if not specified
    end

    % Define time variables from metadata
    tmeta.times = sortrows(tmeta.times); % Sort times if they aren't in order
    % *** NEW: Ensure tmeta.times is numeric (addresses previous type error) ***
    if iscell(tmeta.times)
        try
            tmeta.times = str2double(tmeta.times);
        catch
            error('tmeta.times is a cell array and cannot be converted to numeric. Check its content.');
        end
    end
    % **************************************************************************
    t0 = tmeta.times(1); % Initial time
    tint = mean(diff(tmeta.times)); % Calculate average time interval
    disp("Time steps are : " + tint);
    tf = tmeta.times(end); % Final time
    t = tmeta.times; % Use actual times from tmeta
    tval = t0:0.1:tf; % Finer time points for sine-fitted values

    % Compute circadian information for the custom gene and cell type
    [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, false, period12, ...
                                          cust_gene, cust_cells, false, norm_str);
    
    % Find the specified gene in the results
    gene_idx = find(strcmp(T1.Genes, cust_gene));
    if isempty(gene_idx)
        error('Specified gene not found in the dataset.');
    end

    if strcmp(norm_str, 'lib_size')
        % This is regular cells pipeline
        X_norm = pkg.norm_libsize(sce.X, 1e4); % Renamed X to X_norm to avoid conflict with table variable
        X_norm = log1p(X_norm);
    else % 'magic_impute'
        % This for cancer cells
        X_norm = sc_impute(sce.X, 'MAGIC');
    end
    sce.X = X_norm; % Assign the normalized data back
    clear X_norm;

    % Subset and normalize cells
    ic0 = find(sce.c_cell_type_tx == cust_cells);
    sce_sub = sce.selectcells(ic0);
    ig = find(sce_sub.g == cust_gene);
    clear ic0;

    % Prepare gene expression data for both plot types (scatter and table for violin)
    batch_time = unique(sce_sub.c_batch_id); 
    
    % Pre-allocate based on total cells for the gene within sce_sub
    ncell_total_sub = length(find(sce.c_cell_type_tx == cust_cells)); % Total relevant cells
    
    % These will be populated for both scatter and for creating the table
    all_expression_data = zeros(ncell_total_sub, 1);
    all_time_points = zeros(ncell_total_sub, 1);
    
    ncell_cummu = 0; % Keep track of cumulative cells for indexing
    
    for it = 1:length(batch_time) % Iterate through unique time points (batches)
        current_batch_id = batch_time(it);
        ics = find(sce_sub.c_batch_id == current_batch_id); % Get indices of cells for this time point
        current_expression_data = full(sce_sub.X(ig, ics)); % Get expression for current batch
        nc_loc = length(ics); % Number of cells at this time point
        
        if ig <= size(sce_sub.X, 1) && nc_loc > 0 % Ensure valid gene index and cells exist
            ibeg = ncell_cummu + 1;
            iend = ncell_cummu + nc_loc;
            
            all_expression_data(ibeg:iend) = current_expression_data;
            all_time_points(ibeg:iend) = t(it); % Use actual time from `t`
        end
        ncell_cummu = ncell_cummu + nc_loc; % Accumulate total cell count
    end

    % Trim pre-allocated arrays to actual size if fewer cells were found
    if ncell_cummu < ncell_total_sub
        all_expression_data = all_expression_data(1:ncell_cummu);
        all_time_points = all_time_points(1:ncell_cummu);
    end

    clear sce_sub; % Clear sce_sub after data extraction

    % Create table for violin plot with numeric TimePoints for XData
    tbl_violin = table(all_time_points, all_expression_data, ...
                       'VariableNames', {'TimePoints', 'Expression'});
    
    % Generate the plot for the specified gene if found
    % Set the current axes to the provided axHandle
    axes(axHandle);

    % Generate sine-fitted values based on circadian parameters
    fval = T1.Amp(gene_idx) * cos(2 * pi * (tval - T1.Acrophase(gene_idx)) / T1.Period(gene_idx)) + T1.Mesor(gene_idx);
    Rzts = table2array(T2(gene_idx, 2:end));

    % Plot sine-fitted values
    plot(tval, fval, 'LineWidth', 2);
    hold on;

    % Plot original expression data (Rzts)
    plot(t, Rzts, 'o-', 'LineWidth', 1.5);

    if print_scdata
        use_violin_plot = true;
        if use_violin_plot
            % Plot violin plots from table
            % Using the table syntax: violinplot(T, DataVar, GroupingVar)
            h = violinplot(tbl_violin, "TimePoints", "Expression"); 
            sc_data_legend_str = 'sc-data'; % Update legend string
        else % use_violin_plot is false, so plot scatter
            % Plot scatter points with transparency to highlight overlaps
            scatter(all_time_points, all_expression_data, 20, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
                'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.12); % Set alpha for transparency
            sc_data_legend_str = 'sc-data (Scatter)'; % Update legend string
        end
    end

    % Plot vertical line for phase
    xline(T1.Acrophase_24(gene_idx), 'LineStyle', '--', 'Color', 'r');

    % Set plot limits and labels
    xlim([t0 tf]);
    % Modify the x-axis ticks to show every 1 hour
    xticks(t0:tint:tf);

    % Titles and legends
    title_str = sprintf(['Gene - %s | p-val F-test: %.3f | p-value corr: %.3f \n Phase: %.3f | NumCells: %d'] ,...
                      T1.Genes{gene_idx}, T1.pvalue_corr, T1.pvalue(gene_idx), ...
                      T1.Acrophase_24(gene_idx), ncell_cummu );
    title(title_str, 'FontSize', 16); % Increase title font size
    xlabel('Time (hrs)', 'FontSize', 14); % Increase x-axis label font size
    ylabel('Expression', 'FontSize', 14); % Increase y-axis label font size
    set(gca, 'FontSize', 12); % Increase axis numbers font size

    if print_scdata
        ylim([min(all_expression_data) max(all_expression_data)]);
        legend({'Sine-fitted expression', 'Mean Expression (Rzts)', ...
            'sc-data', 'Phase'}, 'Location', 'northwest');
    else
        legend({'Sine-fitted expression', 'Mean Expression (Rzts)', ...
            'Phase'}, 'Location', 'northwest');
    end
    hold off;
end
