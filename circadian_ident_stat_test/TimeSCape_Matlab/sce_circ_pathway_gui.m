function sce_circ_pathway_gui(sce, tmeta)
%SCE_CIRC_PATHWAY_GUI  Interactive pathway circadian analysis window.
%   sce_circ_pathway_gui(sce, tmeta) opens a window to: load a confident
%   circadian gene list (CSV), bin genes by acrophase (user width), enrich each
%   ZT bin via Enrichr, pick a pathway, and plot its AUCell circadian rhythm
%   (violin + fitted cosine) with a selectable color scheme.
%
%   Requires the TimeSCape MATLAB functions on the path and internet access.

    % ---- shared state ----
    bins     = struct([]);
    pw_map   = zeros(0,2);        % rows -> [bin_index, term_index]

    % ---- figure ----
    fig = figure('Name','TimeSCape — Pathway Circadian Analysis', ...
                 'NumberTitle','off','Color',[1 1 1], ...
                 'Position',[100 100 1060 620]);
    P = @(x,y,w,h) [x y w h];
    lab = @(x,y,w,s) uicontrol('Parent',fig,'Style','text','Position',P(x,y,w,18), ...
              'String',s,'HorizontalAlignment','left','BackgroundColor',[1 1 1],'FontSize',9);

    % ---- controls (left column) ----
    lab(20,582,120,'Confident CSV:');
    hCsv    = uicontrol('Parent',fig,'Style','edit','Position',P(20,560,225,22), ...
                        'HorizontalAlignment','left','String','');
    uicontrol('Parent',fig,'Style','pushbutton','Position',P(250,560,70,24), ...
              'String','Browse','Callback',@onBrowse);

    lab(20,532,120,'Cell type:');
    ctList = cellstr(string(unique(sce.c_cell_type_tx)));
    hCell  = uicontrol('Parent',fig,'Style','popupmenu','Position',P(110,532,210,24), ...
                       'String',ctList);

    lab(20,502,120,'Bin width (h):');
    hBin   = uicontrol('Parent',fig,'Style','edit','Position',P(110,502,60,22),'String','3');

    lab(20,472,120,'Enrichr library:');
    hLib   = uicontrol('Parent',fig,'Style','popupmenu','Position',P(20,448,300,24), ...
                       'String',default_libs());

    uicontrol('Parent',fig,'Style','pushbutton','Position',P(20,410,185,30), ...
              'String','Run Enrichment','FontWeight','bold', ...
              'BackgroundColor',[0.20 0.45 0.70],'ForegroundColor',[1 1 1], ...
              'Callback',@onRunEnrichment);
    hStatus = uicontrol('Parent',fig,'Style','text','Position',P(212,412,108,24), ...
                        'String','','HorizontalAlignment','left','BackgroundColor',[1 1 1]);

    lab(20,380,120,'Pathway (ZT window | term):');
    hPath  = uicontrol('Parent',fig,'Style','popupmenu','Position',P(20,356,300,24), ...
                       'String',{'(run enrichment first)'});

    lab(20,326,120,'Color scheme:');
    hCmap  = uicontrol('Parent',fig,'Style','popupmenu','Position',P(110,326,120,24), ...
                       'String',{'parula','turbo','cool','spring','hot','jet'});
    hViolin= uicontrol('Parent',fig,'Style','checkbox','Position',P(20,296,160,22), ...
                       'String','Violin','Value',1,'BackgroundColor',[1 1 1]);

    uicontrol('Parent',fig,'Style','pushbutton','Position',P(20,258,300,32), ...
              'String','Plot Pathway Rhythm','FontWeight','bold', ...
              'BackgroundColor',[0.15 0.55 0.30],'ForegroundColor',[1 1 1], ...
              'Callback',@onPlot);

    uicontrol('Parent',fig,'Style','pushbutton','Position',P(20,224,300,24), ...
              'String','Refresh full library list (Enrichr)','Callback',@onRefreshLibs);

    % ---- plot axes (right) ----
    hAx = axes('Parent',fig,'Units','pixels','Position',P(370,60,660,520));
    title(hAx,'Run enrichment, choose a pathway, then Plot.','Interpreter','none');

    % ================= callbacks =================
    function onBrowse(~,~)
        [f,p] = uigetfile({'*.csv','CSV files'}, 'Select *_circadian_analysis_confident.csv');
        if isequal(f,0); return; end
        set(hCsv,'String',fullfile(p,f));
    end

    function lib = current_lib()
        items = get(hLib,'String'); v = get(hLib,'Value'); lib = items{v};
        if startsWith(lib,'★') || startsWith(lib,'─')
            error('Pick a library, not a header row.');
        end
    end

    function onRunEnrichment(~,~)
        try
            csv = strtrim(get(hCsv,'String'));
            if isempty(csv) || ~isfile(csv); errordlg('Choose a valid CSV first.'); return; end
            T1c = readtable(csv);
            bw  = str2double(get(hBin,'String')); if isnan(bw)||bw<=0; bw=3; end
            lib = current_lib();
            set(hStatus,'String','binning...'); drawnow;
            bins = sce_circ_bin_genes(T1c, bw);

            labels = {}; pw_map = zeros(0,2);
            for i = 1:numel(bins)
                if numel(bins(i).genes) < 3; bins(i).enrich = []; continue; end
                set(hStatus,'String',sprintf('enrich %s', bins(i).label)); drawnow;
                E = sce_circ_enrichr(bins(i).genes, lib);
                bins(i).enrich = E;
                for k = 1:height(E)
                    labels{end+1} = sprintf('%s | %s  (adjP=%.1e)', ...
                        bins(i).label, E.Term(k), E.AdjPvalue(k)); %#ok<AGROW>
                    pw_map(end+1,:) = [i k]; %#ok<AGROW>
                end
            end
            if isempty(labels); labels = {'(no enriched pathways)'}; end
            set(hPath,'String',labels,'Value',1);
            set(hStatus,'String',sprintf('%d pathways', size(pw_map,1)));
        catch ME
            set(hStatus,'String','error'); errordlg(ME.message,'Enrichment error');
        end
    end

    function onPlot(~,~)
        try
            if isempty(pw_map); errordlg('Run enrichment and pick a pathway first.'); return; end
            sel = get(hPath,'Value');
            i = pw_map(sel,1); k = pw_map(sel,2);
            pw   = bins(i).enrich.Genes{k};
            name = sprintf('%s | %s', bins(i).label, bins(i).enrich.Term(k));
            cts  = get(hCell,'String'); ct = cts{get(hCell,'Value')};
            cms  = get(hCmap,'String'); cmap = cms{get(hCmap,'Value')};
            uv   = get(hViolin,'Value') == 1;
            sce_circ_pathway_plot(sce, tmeta, pw, name, ct, false, hAx, cmap, uv);
        catch ME
            errordlg(ME.message,'Plot error');
        end
    end

    function onRefreshLibs(~,~)
        try
            set(hStatus,'String','fetching libs...'); drawnow;
            menu = sce_circ_enrichr_menu();
            set(hLib,'String',menu,'Value',2);   % row 2 = first suggested
            set(hStatus,'String','libs loaded');
        catch ME
            set(hStatus,'String','lib error'); errordlg(ME.message,'Library list error');
        end
    end
end

function libs = default_libs()
% shown immediately without a network call; use "Refresh" for the full list
    libs = {'GO_Biological_Process_2026','GO_Molecular_Function_2026', ...
            'GO_Cellular_Component_2026','KEGG_2019_Mouse', ...
            'Reactome_Pathways_2024','WikiPathways_2024_Mouse','MSigDB_Hallmark_2020'};
end
