function TimeSCape_GUI
% TIMESCAPE_GUI  –  TimeSCape Circadian Analysis (v0.2)
%
% WORKFLOW FOR BIOLOGISTS:
%   1. The SCE file is selected at startup.
%   2. Click "Define / Edit ZT Times" to map your batch IDs to hours.
%   3. Choose cell type, normalization, period.
%   4. Click "Run Analysis" → produces ALL-genes + confident-only files.
%   5. Click "Generate Heatmap" to visualise the confident genes.
%   6. Use the Gene Explorer panel to plot individual genes.
%
% OUTPUT FILES (written to current MATLAB folder):
%   *_circadian_analysis_all.csv          ← every gene tested
%   *_circadian_analysis_confident.csv    ← F-test AND corr p < 0.05
%   *_circadian_ZTs_mean.csv              ← raw per-ZT means (all)
%   *_circadian_ZTs_mean_normalized.csv   ← ZT00-normalised (all)
%   *_circadian_ZTs_mean_confident.csv    ← per-ZT means (confident)
%   *_circadian_ZTs_mean_normalized_confident.csv

    % ── Style constants ────────────────────────────────────────────────────
    BG       = [1.00 1.00 1.00];
    PAN_BG   = [0.95 0.96 0.98];
    ACCENT   = [0.15 0.32 0.62];
    BTN_FG   = [1.00 1.00 1.00];
    BTN_RUN  = [0.10 0.48 0.24];   % green  – main action
    BTN_PLOT = [0.22 0.44 0.69];   % blue   – plotting
    BTN_WARN = [0.62 0.18 0.18];   % red    – destructive / close
    FNT      = 'Helvetica';
    FSZ      = 11;

    % ── Main figure ────────────────────────────────────────────────────────
    hFig = figure( ...
        'Name',       'TimeSCape  –  Circadian Analysis  (v0.2)', ...
        'Position',   [60, 40, 1180, 760], ...
        'Color',       BG, ...
        'MenuBar',    'none', ...
        'ToolBar',    'none', ...
        'NumberTitle','off', ...
        'DefaultUicontrolFontName', FNT, ...
        'DefaultUicontrolFontSize', FSZ, ...
        'DefaultUicontrolBackgroundColor', BG);

    % ── Load SCE ──────────────────────────────────────────────────────────
    sce     = load_sce_data();
    guiData = struct('tmeta', [], 'outdir', '', 'dark', false);   % shared state

    % ──────────────────────────────────────────────────────────────────────
    %   LEFT CONTROL PANEL
    % ──────────────────────────────────────────────────────────────────────
    cpan = uipanel('Parent', hFig, ...
        'Position',        [0.005, 0.005, 0.295, 0.990], ...
        'BackgroundColor',  PAN_BG, ...
        'BorderType',      'line', ...
        'HighlightColor',  [0.70 0.72 0.78], ...
        'Title',           '  Settings & Controls  ', ...
        'FontName',  FNT,  'FontSize', FSZ+1, ...
        'FontWeight','bold','ForegroundColor', ACCENT);

    % Pixel-coordinate helper  (panel ≈ 345 × 745 px; Units = pixels)
    pn = @(x,y,w,h) [x, y, w, h];

    % ── ① ZT metadata ─────────────────────────────────────────────────────
    section_label(cpan, pn(8,700,320,18), '① Define ZT Time Metadata', ACCENT, PAN_BG);
    mk_btn(cpan, pn(12,670,200,30), 'Define / Edit ZT Times', @defineTmetaCallback, ACCENT, BTN_FG);
    hTmetaStatus = uicontrol('Parent', cpan, 'Style', 'text', ...
        'Position', pn(8,648,320,20), ...
        'String',   '⚠  Tmeta not defined yet', ...
        'ForegroundColor', [0.7 0.3 0], ...
        'BackgroundColor', PAN_BG, ...
        'HorizontalAlignment', 'left', 'FontSize', FSZ);

    % ── ② Analysis settings ────────────────────────────────────────────────
    section_label(cpan, pn(8,618,320,18), '② Analysis Settings', ACCENT, PAN_BG);

    uicontrol('Parent',cpan,'Style','text','Position',pn(10,592,120,20),...
        'String','Cell Type:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    cell_types = unique(sce.c_cell_type_tx);
    hCells = uicontrol('Parent',cpan,'Style','popupmenu',...
        'Position',pn(135,594,195,22),'String',cell_types,'BackgroundColor','white');

    uicontrol('Parent',cpan,'Style','text','Position',pn(10,562,120,20),...
        'String','Normalization:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    hNorm = uicontrol('Parent',cpan,'Style','popupmenu',...
        'Position',pn(135,564,195,22),...
        'String',{'lib_size','magic_impute'},'BackgroundColor','white');

    hPeriod12 = uicontrol('Parent',cpan,'Style','checkbox',...
        'Position',pn(12,534,310,22),...
        'String','Use 12-hr period  (default: 24-hr period)',...
        'BackgroundColor',PAN_BG);

    % ── ③ Run Analysis ─────────────────────────────────────────────────────
    section_label(cpan, pn(8,504,320,18), '③ Run Analysis  →  writes ALL + Confident files', ACCENT, PAN_BG);

    % Brief output description for biologists
    uicontrol('Parent',cpan,'Style','text','Position',pn(12,442,318,58),...
        'String', sprintf(['Produces two result sets:\n' ...
                           '• ALL genes  (every gene tested)\n' ...
                           '• Confident  (F-test AND corr p<0.05)']),...
        'BackgroundColor',PAN_BG,'ForegroundColor',[0.12 0.38 0.12],...
        'HorizontalAlignment','left','FontSize',FSZ);

    mk_btn(cpan, pn(12,408,200,32), '▶  Run Analysis', @runAnalysisCallback, BTN_RUN, BTN_FG);
    mk_btn(cpan, pn(218,408,110,32), '▶▶  All Types', @runAllCallback, [0.06 0.34 0.18], BTN_FG);
    hRunStatus = uicontrol('Parent',cpan,'Style','text',...
        'Position',pn(8,378,320,26),...
        'String','','BackgroundColor',PAN_BG,...
        'HorizontalAlignment','left','FontSize',FSZ,...
        'ForegroundColor',[0.10 0.10 0.10]);

    % ── ④ Heatmap ─────────────────────────────────────────────────────────
    section_label(cpan, pn(8,350,320,18), '④ Generate Heatmap', ACCENT, PAN_BG);

    hStrictFilter = uicontrol('Parent',cpan,'Style','checkbox',...
        'Position',pn(12,324,310,22),...
        'String','Show confident genes only  (recommended)',...
        'Value',1,'BackgroundColor',PAN_BG);
    hClassicCirc  = uicontrol('Parent',cpan,'Style','checkbox',...
        'Position',pn(12,298,280,22),...
        'String','Restrict to core circadian gene set',...
        'Value',0,'BackgroundColor',PAN_BG);

    uicontrol('Parent',cpan,'Style','text','Position',pn(12,272,110,20),...
        'String','Custom label:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    hCustomName = uicontrol('Parent',cpan,'Style','edit',...
        'Position',pn(128,270,200,24),'String','','BackgroundColor','white');

    mk_btn(cpan, pn(12,234,200,30), 'Generate Heatmap', @heatmapCallback, BTN_PLOT, BTN_FG);

    % ── ⑤ Gene explorer ───────────────────────────────────────────────────
    section_label(cpan, pn(8,204,320,18), '⑤ Gene Explorer', ACCENT, PAN_BG);

    plot_types = {'Confident genes','Non-confident genes','Classic circadian genes'};
    uicontrol('Parent',cpan,'Style','text','Position',pn(10,178,110,20),...
        'String','Batch plots:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    hPlotType = uicontrol('Parent',cpan,'Style','popupmenu',...
        'Position',pn(125,180,205,22),'String',plot_types,'BackgroundColor','white');
    mk_btn(cpan, pn(12,146,220,30), 'Save Gene Plots (Batch)', @batchPlotsCallback, BTN_PLOT, BTN_FG);

    uicontrol('Parent',cpan,'Style','text','Position',pn(10,114,130,20),...
        'String','Single gene:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    genes_sorted = sort(sce.g);
    hGeneDropdown = uicontrol('Parent',cpan,'Style','popupmenu',...
        'Position',pn(145,116,185,22),'String',genes_sorted,'BackgroundColor','white');
    uicontrol('Parent',cpan,'Style','text','Position',pn(10,86,90,20),...
        'String','Or type:','BackgroundColor',PAN_BG,'HorizontalAlignment','left');
    hGeneEdit = uicontrol('Parent',cpan,'Style','edit',...
        'Position',pn(105,86,175,24),'String','','BackgroundColor','white');

    hPrintSCdata = uicontrol('Parent',cpan,'Style','checkbox',...
        'Position',pn(12,60,240,22),'String','Overlay single-cell data',...
        'BackgroundColor',PAN_BG);

    % Violin / scatter toggle (only active when overlay is checked)
    uicontrol('Parent',cpan,'Style','text','Position',pn(12,38,80,18),...
        'String','SC style:','BackgroundColor',PAN_BG,'HorizontalAlignment','left',...
        'FontSize',FSZ);
    hViolinBtn = uicontrol('Parent',cpan,'Style','radiobutton',...
        'Position',pn(95,38,82,18),'String','Violin',...
        'BackgroundColor',PAN_BG,'Value',1,'FontSize',FSZ);
    hScatterBtn = uicontrol('Parent',cpan,'Style','radiobutton',...
        'Position',pn(183,38,90,18),'String','Dots',...
        'BackgroundColor',PAN_BG,'Value',0,'FontSize',FSZ);

    % Keep the two radio buttons mutually exclusive
    set(hViolinBtn,  'Callback', @(~,~) set(hScatterBtn,'Value',0));
    set(hScatterBtn, 'Callback', @(~,~) set(hViolinBtn, 'Value',0));

    mk_btn(cpan, pn(12,10,220,26), 'Plot Single Gene', @plotGeneCallback, BTN_PLOT, BTN_FG);

    % ──────────────────────────────────────────────────────────────────────
    %   RIGHT PLOT AREA
    % ──────────────────────────────────────────────────────────────────────
    hPlotAxes = axes('Parent', hFig, ...
        'Position',  [0.315, 0.08, 0.675, 0.87], ...
        'Color',     'white', ...
        'XColor',    [0.2 0.2 0.2], ...
        'YColor',    [0.2 0.2 0.2], ...
        'Box',       'on', ...
        'FontName',  FNT, 'FontSize', 12);
    xlabel(hPlotAxes, 'Zeitgeber Time (hrs)');
    ylabel(hPlotAxes, 'Expression (log-normalised)');
    title(hPlotAxes,  '← Use the Gene Explorer panel to plot a single gene');

    % ── Save figure button (above plot area, right side) ─────────────────
    uicontrol('Parent', hFig, ...
        'Style',           'pushbutton', ...
        'Position',        [920, 730, 130, 24], ...
        'String',          '💾  Save Figure…', ...
        'FontName',        FNT, 'FontSize', 9, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.22 0.44 0.69], ...
        'ForegroundColor', [1 1 1], ...
        'Callback',        @saveFigCallback);

    % ── Dark / Light theme toggle (top-right of figure) ───────────────────
    hThemeBtn = uicontrol('Parent', hFig, ...
        'Style',           'pushbutton', ...
        'Position',        [1058, 730, 114, 24], ...
        'String',          '🌙  Dark Mode', ...
        'FontName',        FNT, 'FontSize', 9, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.28 0.28 0.32], ...
        'ForegroundColor', [0.92 0.92 0.92], ...
        'Callback',        @themeCallback);

    % ====================================================================
    %   CALLBACK FUNCTIONS
    % ====================================================================

    % ── ① Define / Edit Tmeta ─────────────────────────────────────────────
    function defineTmetaCallback(~, ~)
        batches = unique(sce.c_batch_id);
        if isstring(batches);     old_labels = cellstr(batches);
        elseif isnumeric(batches);old_labels = cellstr(num2str(batches));
        else;                      old_labels = batches; end

        nb = numel(old_labels);
        cell_counts = arrayfun(@(i) sum(strcmp(sce.c_batch_id, old_labels{i})), ...
                               1:nb, 'UniformOutput', true);

        % Pre-fill times if batches already look like 'ZTnn'
        prefilled = zeros(nb,1);
        for i = 1:nb
            lbl = old_labels{i};
            if startsWith(lbl,'ZT') && strlength(lbl)==4
                v = str2double(lbl(3:4));
                if ~isnan(v); prefilled(i) = v; end
            end
        end

        init           = cell(nb,3);
        init(:,1)      = old_labels;
        init(:,2)      = num2cell(prefilled);
        init(:,3)      = num2cell(cell_counts(:));

        tfig = figure('Name','Define ZT Tmeta','Position',[180,180,590,450],...
                      'Color',BG,'MenuBar','none','ToolBar','none',...
                      'NumberTitle','off');

        uicontrol('Parent',tfig,'Style','text',...
            'Position',[15,400,560,30],...
            'String',['Set the ZT hour for each batch. ' ...
                      'Enter -1 to exclude a batch completely.'],...
            'BackgroundColor',BG,'HorizontalAlignment','left',...
            'ForegroundColor',[0.10 0.10 0.10]);

        tbl = uitable('Parent',tfig,'Position',[15,80,560,310],...
            'Data',init,...
            'ColumnName',{'Original Batch ID','ZT Hour (numeric)','# Cells'},...
            'ColumnEditable',[false,true,false],...
            'ColumnWidth',{200,150,100},...
            'BackgroundColor',[1 1 1; 0.97 0.97 1]);

        mk_btn(tfig,[15,20,130,40],'Save & Apply',@saveTmeta,ACCENT,BTN_FG);
        mk_btn(tfig,[160,20,130,40],'Save SCE .mat',@saveNewSCE,[0.4 0.4 0.4],BTN_FG);
        mk_btn(tfig,[305,20,110,40],'Close',@(~,~)close(tfig),BTN_WARN,BTN_FG);

        function saveTmeta(~,~)
            d          = get(tbl,'Data');
            old_lbls   = d(:,1);
            times      = cell2mat(d(:,2));
            cnts       = cell2mat(d(:,3));

            % Build new ZT labels
            new_lbls = cell(numel(times),1);
            for ii = 1:numel(times)
                if times(ii) >= 0; new_lbls{ii} = sprintf('ZT%02d', times(ii));
                else;              new_lbls{ii} = '-1'; end
            end

            % Build tmeta with canonical column names expected by pipeline
            old_labels = old_lbls;
            new_labels = new_lbls;
            % Update these lines in TimeSCape_GUI -> saveTmeta
            guiData.tmeta = sortrows( ...
                table(old_labels, new_labels, times, ...
                      'VariableNames', {'old_labels','new_labels','ZT_times'}), ...
                'ZT_times');

            % Remove excluded batches from working sce
            excluded = old_labels(strcmp(new_labels,'-1'));
            for ib = 1:length(excluded)
                sce.c_batch_id(strcmp(sce.c_batch_id, excluded{ib})) = "-1";
            end
            keep      = ~strcmp(sce.c_batch_id,"-1");
            sce_new   = SingleCellExperiment(sce.X(:,keep), sce.g);
            sce_new.c_cell_type_tx = sce.c_cell_type_tx(keep);
            sce_new.c_cell_id      = sce.c_cell_id(keep);
            sce_new.c_batch_id     = sce.c_batch_id(keep);

            % Merge batches sharing the same new label
            [uniq,~,li] = unique(guiData.tmeta.new_labels);
            for ii = 1:length(uniq)
                m = ismember(sce_new.c_batch_id, guiData.tmeta.old_labels(li==ii));
                sce_new.c_batch_id(m) = uniq{ii};
            end

            clear sce; sce = sce_new; clear sce_new;
            assignin('base','sce',sce);

            set(hTmetaStatus,'String', ...
                sprintf('✓  Tmeta set: %d time points', sum(times >= 0)), ...
                'ForegroundColor',[0.1 0.5 0.1]);
            disp('Tmeta saved and SCE updated.');
        end

        function saveNewSCE(~,~)
            [file,path] = uiputfile('*.mat','Save Modified SCE As');
            if ~isequal(file,0)
                save(fullfile(path,file),'sce');
                disp(['SCE saved → ', fullfile(path,file)]);
            end
        end
    end

    % ── ③ Run Analysis ─────────────────────────────────────────────────────
    function runAnalysisCallback(~,~)
        if ~check_tmeta(); return; end
        celltype = get_sel(hCells);
        period12 = logical(hPeriod12.Value);
        norm_str = get_sel(hNorm);

        set(hRunStatus,'String','⏳  Running analysis…','ForegroundColor',[0.4 0.2 0]);
        drawnow;
        try
            [T1, T2] = sce_circ_phase_estimation_stattest(sce, guiData.tmeta, ...
                           true, period12, [], celltype, false, norm_str);
            n_all  = height(T1);
            n_conf = sum((T1.pvalue < 0.05) & (T1.pvalue_corr < 0.05));
            % Store output directory so other callbacks can find the files
            outdir_name        = regexprep(strtrim(char(celltype)), '[^\w]', '_');
            guiData.outdir     = fullfile(pwd, outdir_name);
            set(hRunStatus,'String', ...
                sprintf('✓  Done: %d genes tested, %d confident  →  %s', ...
                        n_all, n_conf, outdir_name), ...
                'ForegroundColor',[0.1 0.5 0.1]);
        catch ME
            set(hRunStatus,'String','✗  Error – see Command Window','ForegroundColor',[0.7 0 0]);
            errordlg(['Analysis error: ', ME.message], 'Error');
        end
    end

    % ── ③b Run ALL cell types ─────────────────────────────────────────────
    function runAllCallback(~,~)
        if ~check_tmeta(); return; end
        period12 = logical(hPeriod12.Value);
        norm_str = get_sel(hNorm);

        all_types = unique(sce.c_cell_type_tx);
        n_types   = numel(all_types);

        set(hRunStatus,'String', ...
            sprintf('⏳  Running all %d cell types…', n_types), ...
            'ForegroundColor',[0.4 0.2 0]);
        drawnow;

        n_all_total  = 0;
        n_conf_total = 0;
        errors       = {};

        for ii = 1:n_types
            ct = char(all_types(ii));
            set(hRunStatus,'String', ...
                sprintf('⏳  [%d/%d]  %s…', ii, n_types, ct), ...
                'ForegroundColor',[0.4 0.2 0]);
            drawnow;
            try
                [T1, ~] = sce_circ_phase_estimation_stattest(sce, guiData.tmeta, ...
                               true, period12, [], ct, false, norm_str);
                n_all_total  = n_all_total  + height(T1);
                n_conf_total = n_conf_total + sum((T1.pvalue < 0.05) & (T1.pvalue_corr < 0.05));
            catch ME
                errors{end+1} = sprintf('%s: %s', ct, ME.message); %#ok<AGROW>
            end
        end

        if isempty(errors)
            set(hRunStatus,'String', ...
                sprintf('✓  All %d types done — %d genes, %d confident', ...
                        n_types, n_all_total, n_conf_total), ...
                'ForegroundColor',[0.1 0.5 0.1]);
        else
            set(hRunStatus,'String', ...
                sprintf('⚠  Done with %d error(s) — see Command Window', numel(errors)), ...
                'ForegroundColor',[0.65 0.3 0]);
            fprintf('\n--- Run-All errors ---\n');
            for e = errors; disp(e{1}); end
        end
        % Point outdir to the last successfully processed cell type directory
        outdir_name    = regexprep(strtrim(char(all_types(end))), '[^\w]', '_');
        guiData.outdir = fullfile(pwd, outdir_name);
    end

    % ── ④ Heatmap ─────────────────────────────────────────────────────────
    function heatmapCallback(~,~)
        if ~check_tmeta(); return; end
        celltype    = get_sel(hCells);
        period12    = logical(hPeriod12.Value);
        strict      = logical(hStrictFilter.Value);
        circ        = logical(hClassicCirc.Value);
        customName  = strtrim(hCustomName.String);
        if period12; pl = 'period_12'; else; pl = 'period_24'; end

        fname_chk = sprintf('%s_%s_circadian_analysis_all.csv', celltype, pl);
        if ~exist(fname_chk,'file')
            warndlg(sprintf(['Analysis file not found:\n%s\n\n' ...
                             'Please run the analysis first (Step ③).'], fname_chk), ...
                    'File Not Found');
            return;
        end
        try
            generateHeatmap_circ_simple(celltype, strict, customName, circ, period12, ...
                                        guiData.outdir, hPlotAxes);
            msgbox(sprintf('Heatmap saved →  %s', guiData.outdir), 'Done');
        catch ME
            errordlg(['Heatmap error: ', ME.message], 'Error');
        end
    end

    % ── ⑤a Batch gene plots ───────────────────────────────────────────────
    function batchPlotsCallback(~,~)
        if ~check_tmeta(); return; end
        cust_cells = get_sel(hCells);
        plot_type  = hPlotType.Value;
        period12   = logical(hPeriod12.Value);
        norm_str   = get_sel(hNorm);
        ws = warning('off','MATLAB:mkdir:DirectoryExists');
        try
            sce_circ_plot(sce, guiData.tmeta, cust_cells, plot_type, period12, norm_str, guiData.outdir);
            msgbox(sprintf('Gene plots saved →  %s', guiData.outdir),'Done');
        catch ME
            errordlg(['Plot error: ', ME.message],'Error');
        end
        warning(ws);
    end

    % ── ⑤b Single gene plot ───────────────────────────────────────────────
    function plotGeneCallback(~,~)
        if ~check_tmeta(); return; end

        typed = strtrim(get(hGeneEdit,'String'));
        if ~isempty(typed); gene = typed;
        else
            opts = get(hGeneDropdown,'String');
            gene = opts{get(hGeneDropdown,'Value')};
        end

        if isempty(gene) || ~ismember(gene, sce.g)
            errordlg(sprintf('Gene "%s" not found in dataset.\nCheck spelling.',gene),...
                     'Invalid Gene');
            return;
        end

        cust_cells    = get_sel(hCells);
        period12      = logical(hPeriod12.Value);
        print_scdata  = logical(hPrintSCdata.Value);
        norm_str      = get_sel(hNorm);
        use_violin    = logical(hViolinBtn.Value);   % true=violin, false=dots

        ws = warning('off','MATLAB:mkdir:DirectoryExists');
        try
            cla(hPlotAxes);
            sce_circ_plot_gene(sce, guiData.tmeta, cust_cells, period12, ...
                               gene, hPlotAxes, print_scdata, norm_str, use_violin, guiData.outdir);
        catch ME
            errordlg(['Plot error: ', ME.message],'Error');
        end
        warning(ws);
    end

    % ====================================================================
    %   SAVE FIGURE
    % ====================================================================

    function saveFigCallback(~,~)
        % Let user choose filename + format
        [fname, fpath] = uiputfile( ...
            {'*.png', 'PNG image  (*.png)'; ...
             '*.svg', 'SVG vector (*.svg)'; ...
             '*.pdf', 'PDF vector (*.pdf)'; ...
             '*.fig', 'MATLAB figure (*.fig)'}, ...
            'Save Current Plot As', ...
            'TimeSCape_plot.png');

        if isequal(fname, 0); return; end   % user cancelled

        fullpath = fullfile(fpath, fname);
        [~, ~, ext] = fileparts(fname);

        % Copy axes (+ colorbar if present) into a clean white export figure
        hExp = figure('Color','white','Visible','off', ...
                      'Position',[0 0 1000 650]);
        newAx = copyobj(hPlotAxes, hExp);
        newAx.Units    = 'normalized';
        newAx.Position = [0.10, 0.11, 0.82, 0.80];

        try
            switch lower(ext)
                case '.fig'
                    saveas(hExp, fullpath);
                case {'.pdf','.svg'}
                    exportgraphics(hExp, fullpath, ...
                        'ContentType','vector', 'BackgroundColor','white');
                otherwise   % .png (default)
                    exportgraphics(hExp, fullpath, ...
                        'Resolution', 200, 'BackgroundColor','white');
            end
            msgbox(sprintf('Saved:\n%s', fullpath), 'Figure Saved', 'help');
        catch ME
            errordlg(['Save failed: ' ME.message], 'Save Error');
        end
        close(hExp);
    end

    % ====================================================================
    %   THEME TOGGLE
    % ====================================================================

    function themeCallback(~,~)
        guiData.dark = ~guiData.dark;
        if guiData.dark
            fig_bg  = [0.12 0.12 0.16];
            pan_bg  = [0.17 0.17 0.22];
            txt_col = [0.88 0.88 0.90];
            acc_col = [0.50 0.72 1.00];
            ax_bg   = [0.10 0.10 0.13];
            ax_tc   = [0.80 0.80 0.82];
            grid_c  = [0.28 0.28 0.32];
            inp_bg  = [0.22 0.22 0.28];
            set(hThemeBtn,'String','☀  Light Mode');
        else
            fig_bg  = [1.00 1.00 1.00];
            pan_bg  = [0.95 0.96 0.98];
            txt_col = [0.10 0.10 0.10];
            acc_col = [0.15 0.32 0.62];
            ax_bg   = [1.00 1.00 1.00];
            ax_tc   = [0.15 0.15 0.15];
            grid_c  = [0.82 0.82 0.82];
            inp_bg  = [1.00 1.00 1.00];
            set(hThemeBtn,'String','🌙  Dark Mode');
        end

        % ── Figure & panel ────────────────────────────────────────────────
        set(hFig,  'Color', fig_bg);
        set(cpan,  'BackgroundColor', pan_bg, 'ForegroundColor', acc_col, ...
                   'HighlightColor', ax_tc*0.6);

        % ── Plot axes ─────────────────────────────────────────────────────
        set(hPlotAxes, 'Color', ax_bg, 'XColor', ax_tc, 'YColor', ax_tc, ...
                       'GridColor', grid_c, 'GridAlpha', 0.7);
        set(get(hPlotAxes,'Title'),  'Color', ax_tc);
        set(get(hPlotAxes,'XLabel'), 'Color', ax_tc);
        set(get(hPlotAxes,'YLabel'), 'Color', ax_tc);

        % ── All uicontrols inside panel ───────────────────────────────────
        kids = findobj(cpan, 'Type','uicontrol');
        for k = 1:numel(kids)
            sty = lower(get(kids(k),'Style'));
            switch sty
                case {'text','checkbox','radiobutton'}
                    set(kids(k), 'BackgroundColor', pan_bg, ...
                                 'ForegroundColor', txt_col);
                case {'popupmenu','edit'}
                    set(kids(k), 'BackgroundColor', inp_bg, ...
                                 'ForegroundColor', txt_col);
            end
        end

        % ── Section labels get accent colour (bold text) ──────────────────
        sects = findobj(cpan,'Type','uicontrol','Style','text','FontWeight','bold');
        for k = 1:numel(sects)
            set(sects(k), 'ForegroundColor', acc_col);
        end

        % ── Keep semantic status text colours ─────────────────────────────
        % (hTmetaStatus / hRunStatus already track their own colour on update)

        drawnow;
    end

    % ====================================================================
    %   UTILITY FUNCTIONS
    % ====================================================================

    function ok = check_tmeta()
        ok = ~isempty(guiData.tmeta);
        if ~ok
            errordlg(['Please define ZT times first (Step ①).'],'Tmeta Required');
        end
    end

    function str = get_sel(hctrl)
        opts = cellstr(hctrl.String);
        str  = opts{hctrl.Value};
    end

    function sce = load_sce_data()
        [fn,fp] = uigetfile('*.mat','Select SCE .mat file');
        if isequal(fn,0); error('No file selected.'); end
        loaded = load(fullfile(fp,fn));
        if isfield(loaded,'sce'); sce = loaded.sce;
        else; error('The .mat file does not contain a variable named "sce".'); end
        fprintf('SCE loaded: %d genes × %d cells\n', size(sce.X,1), size(sce.X,2));
    end

    function section_label(parent, pos, txt, col, bg)
        uicontrol('Parent',parent,'Style','text','Position',pos,...
            'String',txt,'FontWeight','bold','ForegroundColor',col,...
            'BackgroundColor',bg,'HorizontalAlignment','left',...
            'FontName',FNT,'FontSize',FSZ);
    end

    function btn = mk_btn(parent, pos, lbl, cb, bg, fg)
        btn = uicontrol('Parent',parent,'Style','pushbutton',...
            'Position',pos,'String',lbl,'Callback',cb,...
            'BackgroundColor',bg,'ForegroundColor',fg,...
            'FontName',FNT,'FontSize',FSZ,'FontWeight','bold');
    end

end  % TimeSCape_GUI
  