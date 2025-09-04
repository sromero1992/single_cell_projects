function circadian_gui
    % Create the main figure window
    hFig = figure('Name', 'SCE-Circadian Analysis GUI', 'Position', [100, 100, 1000, 700], ... % Increased figure height
        'DefaultUicontrolFontSize', 12, 'DefaultTextFontSize', 12);

    % Load the sce.mat file
    sce = load_sce_data();

    % Initialize GUI data for tmeta
    guiData = struct('tmeta', [], 'cust_cells', '', 'plot_type', 1, 'gene', '');

    % Create button to define tmeta interactively
    uicontrol('Style', 'pushbutton', 'Position', [50, 600, 150, 30], ... % Adjusted position
              'String', 'Define Tmeta', 'Callback', @defineTmetaAndPlotCallback);

    % Normalization type
    norm_types = {'lib_size', 'magic_impute'};
    % Create UI controls
    uicontrol('Style', 'text', 'Position', [45, 550, 110, 20], 'String', 'Normalization:'); % Adjusted position
    hNorm = uicontrol('Style', 'popupmenu', 'Position', [160, 550, 120, 20], ... % Adjusted position
                       'String', norm_types);

    % Create dropdown for selecting cell type
    cell_types = unique(sce.c_cell_type_tx); % Ensure this matches the field in the 'sce' structure
    % Create UI controls
    uicontrol('Style', 'text', 'Position', [35, 500, 100, 20], 'String', 'Cell Type:'); % Adjusted position
    hCells = uicontrol('Style', 'popupmenu', 'Position', [120, 500, 150, 20], ... % Adjusted position
                       'String', cell_types, 'Value', 1); % Set default selection to the first option

    % Define your custom plot types
    plot_types = {'Plot confident genes', 'Plot non-condifent genes', ...
                  'Plot circadian only', 'Plot syncronized genes'};
    % Create dropdown for selecting plot_type
    uicontrol('Style', 'text', 'Position', [35, 450, 100, 20], 'String', 'Plot Type:'); % Adjusted position
    hPlotType = uicontrol('Style', 'popupmenu', 'Position', [120, 450, 170, 20], ... % Adjusted position
                          'String', plot_types);

    % Create checkbox for selecting period12
    hPeriod12 = uicontrol('Style', 'checkbox', 'Position', [55, 410, 210, 20], ... % Adjusted position
                          'String', 'Period 12 (otherwise 24)');

    % Heatmap Options
    uicontrol('Style', 'text', 'Position', [15, 380, 200, 20], 'String', 'Heatmap Options:'); % Added section header
    hStrictFilter = uicontrol('Style', 'checkbox', 'Position', [55, 350, 280, 20], ... % Added Strict filter checkbox
                              'String', 'Only confident genes', 'Value', 1); % Default to strict
    hClassicCirc = uicontrol('Style', 'checkbox', 'Position', [55, 320, 250, 20], ... % Added Classic circ genes checkbox
                             'String', 'Circadian genes', 'Value', 0); % Default to not filter
    uicontrol('Style', 'text', 'Position', [230, 350, 120, 20], 'String', 'Heatmap Name:'); % Added Custom name label
    hCustomName = uicontrol('Style', 'edit', 'Position', [230, 320, 120, 25], ... % Added Custom name text field
                            'String', ''); % Default empty

    % Create button for Analysis (now also generates heatmap)
    uicontrol('Style', 'pushbutton', 'Position', [50, 270, 170, 30], ... % Adjusted position
              'String', 'Analyze & Heatmap', 'Callback', @analyzeAndHeatmapCallback); % Renamed callback and string

    % Create button for Plot All
    uicontrol('Style', 'pushbutton', 'Position', [50, 230, 170, 30], ... % Adjusted position
              'String', 'Plot Genes & Analyze', 'Callback', @plotAllGenesCallback);

    % Create text box for gene input
    genes = sce.g; % Assuming sce.g is a cell array or similar
    uicontrol('Style', 'text', 'Position', [45, 170, 100, 20], 'String', 'Single gene plot:'); % Adjusted position
    hGene = uicontrol('Style', 'popupmenu', 'Position', [140, 170, 150, 20], ... % Adjusted position
                       'String', genes); % Dropdown menu for genes

    % Create button for Plot Gene
    uicontrol('Style', 'pushbutton', 'Position', [45, 120, 150, 30], ... % Adjusted position
              'String', 'Plot single Gene', 'Callback', @plotGeneCallback);

    % Create checkbox for printing sc-data
    hPrintSCdata = uicontrol('Style', 'checkbox', 'Position', [200, 125, 120, 20], ... % Adjusted position
                          'String', 'Print sc-data');

    % Create axes for plotting
    hPlotAxes = axes('Parent', hFig, 'Position', [0.4, 0.1, 0.55, 0.8]);
    % Set default font size for axes labels and title
    set(hPlotAxes, 'FontSize', 12);

    function defineTmetaAndPlotCallback(~, ~)
        % Get unique batches from sce.c_batch_id
        batches = unique(sce.c_batch_id);

        % Check the format of batches and convert accordingly
        if isnumeric(batches)
            old_labels = cellstr(num2str(batches)); % Convert numeric batches to strings
        elseif iscellstr(batches)
            old_labels = batches; % Already in cell array of strings
        elseif isstring(batches)
            old_labels = cellstr(batches); % Convert string array to cell array of character vectors
        else
            error('Unexpected format for batches.');
        end

        % Calculate the number of cells for each batch
        num_batches = numel(old_labels);
        cell_counts = zeros(num_batches, 1); % Initialize cell counts array
        for i = 1:num_batches
            cell_counts(i) = sum(strcmp(sce.c_batch_id, old_labels{i}));
        end

        % Create a new figure for tmeta definition
        tmetaFig = figure('Name', 'Define tmeta', 'Position', [150, 150, 550, 400]);

        % Initialize table data with old_labels, times, and cell_counts
        initial_data = cell(num_batches, 3);
        initial_data(:, 1) = old_labels; % Set old_labels
        initial_data(:, 3) = num2cell(cell_counts); % Set cell counts for each batch
        
        % --- NEW LOGIC TO PRE-FILL ZT TIMES ---
        prefilled_times = zeros(num_batches, 1); % Default to 0
        for i = 1:num_batches
            current_label = old_labels{i};
            % Check if the label starts with 'ZT' and has 4 characters (e.g., 'ZT00', 'ZT03')
            if startsWith(current_label, 'ZT') && strlength(current_label) == 4
                zt_str = current_label(3:4); % Extract the numeric part (e.g., '00', '03')
                zt_val = str2double(zt_str);
                if ~isnan(zt_val) % Ensure it's a valid number
                    prefilled_times(i) = zt_val;
                end
            end
        end
        initial_data(:, 2) = num2cell(prefilled_times);


        % Create a table UI to input new labels and times, with cell counts as non-editable
        tmetaTable = uitable('Parent', tmetaFig, 'Position', [25, 75, 500, 250], ...
                             'Data', initial_data, ...
                             'ColumnName', {'Old Labels', 'Times', 'Cell Count'}, ...
                             'ColumnEditable', [false, true, false]);  % Cell count column is non-editable

        % Button to save tmeta changes and update SCE
        uicontrol('Style', 'pushbutton', 'Position', [60, 25, 100, 30], ...
                  'String', 'Save tmeta', 'Callback', @saveTmeta);

        % Button to save the modified SCE object
        uicontrol('Style', 'pushbutton', 'Position', [190, 25, 100, 30], ...
                  'String', 'Save SCE', 'Callback', @saveNewSCE);

        % Button to close the tmeta definition window
        uicontrol('Style', 'pushbutton', 'Position', [320, 25, 100, 30], ...
                  'String', 'Close', 'Callback', @(~, ~) close(tmetaFig));

        % Callback to save tmeta and update the SCE object
        function saveTmeta(~, ~)
            % Retrieve the updated table data
            tableData = get(tmetaTable, 'Data');
            old_labels = tableData(:, 1); % This should be cell array of strings
            times = cell2mat(tableData(:, 2)); % Convert times from cell array to numeric
            cell_counts = cell2mat(tableData(:, 3)); % Cell counts remain unchanged

            % Convert times to ZT labels
            new_labels = cell(numel(times), 1);
            for i = 1:numel(times)
                if times(i) >= 0
                    new_labels{i} = sprintf('ZT%02d', times(i));
                else
                    new_labels{i} = sprintf('%d', times(i));
                end
            end

            % Convert to table and sort times (if not sorted)
            tbl = table(old_labels, new_labels, times);
            tbl = sortrows(tbl, 'times');
            guiData.tmeta = tbl; % Update GUI data

            % Process to remove selected batches
            rm_batch = ismember(guiData.tmeta.new_labels, '-1');
            rm_batch = guiData.tmeta.old_labels(rm_batch);

            for ib = 1:length(rm_batch)
                % Use pre-selected batches to remove
                idx_rm = find(strcmp(sce.c_batch_id, rm_batch{ib}));
                % Rename batch
                sce.c_batch_id(idx_rm) = "-1";
            end

            % Subsetting sce dynamically
            idx = ~strcmp(sce.c_batch_id, "-1");
            sce_new = SingleCellExperiment(sce.X(:, idx), sce.g);
            sce_new.c_cell_type_tx = sce.c_cell_type_tx(idx);
            sce_new.c_cell_id = sce.c_cell_id(idx);
            sce_new.c_batch_id = sce.c_batch_id(idx);

            % Merge the new labels that are the same
            [unique_labels, ~, label_idx] = unique(guiData.tmeta.new_labels);
            for i = 1:length(unique_labels)
                current_label = unique_labels{i};
                % Find indices in the new SCE that match this new label
                merge_idx = ismember(sce_new.c_batch_id, guiData.tmeta.old_labels(label_idx == i));
                % Update batch ID with the merged label
                sce_new.c_batch_id(merge_idx) = current_label;
            end

            % Update the tmeta with the new old_labels (now as new_labels)
            unique_new_labels = unique(guiData.tmeta.new_labels);  % Get unique new labels after merging
            updated_times = zeros(length(unique_new_labels), 1); % Initialize updated times
            updated_counts = zeros(length(unique_new_labels), 1); % Initialize updated cell counts

            for i = 1:length(unique_new_labels)
                % Get index of the old label corresponding to the unique new label
                idx_old_label = find(strcmp(guiData.tmeta.new_labels, unique_new_labels{i}), 1);
                updated_times(i) = guiData.tmeta.times(idx_old_label); % Assign corresponding times
                updated_counts(i) = cell_counts(idx_old_label); % Keep cell counts consistent
            end

            updatedData = [unique_new_labels, num2cell(updated_times), num2cell(updated_counts)];
            set(tmetaTable, 'Data', updatedData);  % Update the data in the table

            % Replace the original sce with the new sce
            clear sce; % Clear the old SCE object to save memory
            sce = sce_new; % Assign the new SCE object to 'sce'
            clear sce_new; % Clear the temporary new SCE object

            % Update the GUI with the new SCE object
            assignin('base', 'sce', sce);

            % Notify user of successful update
            disp('tmeta saved and SCE updated.');
        end

        % Callback to save the new SCE object
        function saveNewSCE(~, ~)
            % Prompt user to select a file to save the new SCE object
            [file, path] = uiputfile('*.mat', 'Save Modified SCE As');
            if isequal(file, 0)
                disp('User canceled the file save.');
            else
                save(fullfile(path, file), 'sce');
                disp(['Modified SCE saved to ', fullfile(path, file)]);
            end
        end
    end

    % Callback function for Plot Gene button
    function plotGeneCallback(~, ~)
        % Retrieve parameters
        cust_cells_options = cellstr(hCells.String); % Use cellstr for robustness
        cust_cells = cust_cells_options{hCells.Value};
        period12 = hPeriod12.Value;
        gene_options = hGene.String; % Get gene options
        gene = gene_options{hGene.Value}; % Retrieve selected gene
        print_scdata = hPrintSCdata.Value;
        norm_str_opts = cellstr(hNorm.String);
        norm_str = norm_str_opts{hNorm.Value};

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end

        % Check if gene is empty or not in the dataset
        if isempty(gene) || ~ismember(gene, sce.g)
            errordlg('Selected gene is invalid or not found in the dataset.');
            return;
        end

        % Create and display waitbar
        %hWaitbar = waitbar(0, 'Plotting Gene...');
        disp("Working on Plot Gene...")
        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');

        try
            % Clear the existing plot in the axes
            cla(hPlotAxes);

            % Call the gene plotting function directly with the gene
            sce_circ_plot_gene(sce, guiData.tmeta, cust_cells, period12, ...
                              gene, hPlotAxes, print_scdata, norm_str);

            % Display completion message
            msgbox('Gene plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurredin plot gene: ', ME.message]);
        end

        % Close the waitbar
        %close(hWaitbar);
        disp("Finished Plot Gene")

        % Restore the warning state
        warning(warningState);
    end
    % Callback function for Plot All button
    function plotAllGenesCallback(~, ~)
        % Retrieve parameters
        cust_cells_options = cellstr(hCells.String); % Use cellstr for robustness
        cust_cells = cust_cells_options{hCells.Value};
        %if contains(cust_cells,"None"); cust_cells = []; end
        plot_type = hPlotType.Value;
        period12 = hPeriod12.Value;
        norm_str_opts = cellstr(hNorm.String);
        norm_str = norm_str_opts{hNorm.Value};

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end
        % Create and display waitbar
        %hWaitbar = waitbar(0, 'Plotting Genes...');
        disp("Working on Plot Genes...")
        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');
        try
            % Call the sce_circ_plot function for all genes
            sce_circ_plot(sce, guiData.tmeta, cust_cells, plot_type, period12, norm_str);

            % Display completion message
            msgbox('All plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end
        % Close the waitbar
        %close(hWaitbar);
        disp("Finished Plot Genes!")
        % Restore the warning state
        warning(warningState);
    end
    % Callback function for Analyze and Heatmap button
    function analyzeAndHeatmapCallback(~, ~)
        % Retrieve parameters
        cust_cells_options = cellstr(hCells.String); % Use cellstr to be robust
        celltype = cust_cells_options{hCells.Value}; % Get selected cell type (renamed to celltype for heatmap function)
        period12 = hPeriod12.Value; % Use this if the analysis function needs it
        norm_str_opts = cellstr(hNorm.String);
        norm_str = norm_str_opts{hNorm.Value};

        % Heatmap specific parameters
        strict = hStrictFilter.Value; % Get state of Strict filter checkbox
        circ = hClassicCirc.Value;   % Get state of Classic circ genes checkbox
        customName = hCustomName.String; % Get text from Custom name field

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before analyzing.');
            return;
        end

        % Create and display waitbar (optional, consider for long analysis)
        %hWaitbar = waitbar(0, 'Analyzing and Generating Heatmap...');
        disp("Working on Analysis and Heatmap...")
        disp(['Cell type selected : ', celltype]);
        disp(['Strict filtering: ', string(strict)]);
        disp(['Classic circ genes filter: ', string(circ)]);
        disp(['Custom name: ', customName]);
        disp(['Normalization: ', norm_str]);

        try
            % --- Analysis Step ---
            % Call the sce_circ_phase_estimation function
            % NOTE: This function is expected to generate the CSV file
            % named like strcat(celltype, "_period_24__macro_circadian_analysis.csv")
            % which is required by generateHeatmap_circ_simple.
            [T1, T2] = sce_circ_phase_estimation_stattest(sce, guiData.tmeta, ...
                                     true, period12, [], celltype, false, norm_str);

            % --- Heatmap Generation Step ---
            % Check if the analysis file was created (optional but good practice)
            analysis_fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
            if exist(analysis_fname, 'file')
                disp('Analysis file found. Generating heatmap...');
                % Call the heatmap generation function
                generateHeatmap_circ_simple(celltype, strict, customName, circ);
                disp('Heatmap generated.');
            else
                 warning('Analysis output file "%s" not found. Heatmap cannot be generated.', analysis_fname);
                 msgbox('Analysis completed, but heatmap file not found.', 'Warning', 'warn');
            end


            % Display analysis completion message
             msgbox('Analysis and Heatmap process completed.');

        catch ME
            % Display an error message
            errordlg(['Error occurred during analysis or heatmap generation: ', ME.message]);
        end

        % Close the waitbar (if used)
        %close(hWaitbar);
        disp("Finished Analysis and Heatmap...")
    end
    % Function to load sce data (customize as needed)
    function sce = load_sce_data()
        [fileName, filePath] = uigetfile('*.mat', 'Select SCE Data File');
        if isequal(fileName, 0)
            error('No file selected.');
        end
        loadedData = load(fullfile(filePath, fileName));
        disp("Successfully loaded sce...");

        % Display the fields of the loaded data
        disp('Fields in loaded data:');
        disp(fieldnames(loadedData));

        % Assuming 'sce' is a field in the loaded data
        if isfield(loadedData, 'sce')
             sce = loadedData.sce;
        else
             error('The loaded .mat file does not contain a variable named ''sce''.');
        end


        % Display fields of the 'sce' structure
        disp('Fields in the sce structure:');
        disp(fieldnames(sce));
    end

   
end