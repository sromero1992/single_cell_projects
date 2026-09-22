function [library, bin_width, norm_str, standardize, useViolin, ok] = tscp_config_dialog(info)
%TSCP_CONFIG_DIALOG  Modal popup: Enrichr database, 24h bin width, score mode.
%   Score mode maps to (norm_str, standardize):
%     1  raw counts -> AUCell            ('none',     false)
%     2  lib_size + log1p -> AUCell      ('lib_size', false)   [default]
%     3  standardized (z-score AUCell)   ('lib_size', true)
%   (AUCell is rank-based, so 1 and 2 give the same score; 3 re-centres it.)
%   Colour is chosen in the main window. Returns ok=false on cancel.
    if nargin < 1; info = []; end
    library = ''; bin_width = 3; norm_str = 'lib_size'; standardize = false; useViolin = true; ok = false;
    libs = {'GO_Biological_Process_2026','GO_Molecular_Function_2026', ...
            'GO_Cellular_Component_2026','KEGG_2019_Mouse', ...
            'Reactome_Pathways_2024','WikiPathways_2024_Mouse','MSigDB_Hallmark_2020'};
    defLib = 1; defBin = '3';
    if ~isempty(info)
        if any(strcmp(libs, info.library)); defLib = find(strcmp(libs, info.library), 1);
        else; libs = [{char(info.library)}, libs]; defLib = 1; end
        defBin = num2str(info.bin_width);
    end
    d = dialog('Position',[500 400 380 250],'Name','Pathway options');
    uicontrol('Parent',d,'Style','text','Position',[20 210 340 18], ...
              'String','Enrichr database:','HorizontalAlignment','left');
    hL = uicontrol('Parent',d,'Style','popupmenu','Position',[20 188 340 24],'String',libs,'Value',defLib);
    uicontrol('Parent',d,'Style','text','Position',[20 156 210 18], ...
              'String','Bin width (h) across the 24 h window:','HorizontalAlignment','left');
    hB = uicontrol('Parent',d,'Style','edit','Position',[236 156 50 22],'String',defBin);
    uicontrol('Parent',d,'Style','text','Position',[20 122 340 18], ...
              'String','Score / normalization:','HorizontalAlignment','left');
    hS = uicontrol('Parent',d,'Style','popupmenu','Position',[20 100 340 24], ...
                   'String',{'raw counts -> AUCell', ...
                             'lib_size + log1p -> AUCell', ...
                             'standardized (z-score AUCell)'}, 'Value',2);
    hV = uicontrol('Parent',d,'Style','checkbox','Position',[20 70 140 22],'String','Violin','Value',1);
    uicontrol('Parent',d,'Style','pushbutton','Position',[80 24 100 32],'String','OK', ...
              'FontWeight','bold','Callback',@(~,~) onok());
    uicontrol('Parent',d,'Style','pushbutton','Position',[200 24 100 32],'String','Cancel', ...
              'Callback',@(~,~) delete(d));
    uiwait(d);
    function onok()
        items = get(hL,'String'); library = items{get(hL,'Value')};
        bin_width = str2double(get(hB,'String'));
        if isnan(bin_width) || bin_width <= 0; bin_width = 3; end
        switch get(hS,'Value')
            case 1; norm_str = 'none';     standardize = false;
            case 2; norm_str = 'lib_size'; standardize = false;
            case 3; norm_str = 'lib_size'; standardize = true;
        end
        useViolin = get(hV,'Value') == 1;
        ok = true; delete(d);
    end
end
