function [menu, suggested, allLibs] = sce_circ_enrichr_menu()
%SCE_CIRC_ENRICHR_MENU  Ordered Enrichr library list for the GUI dropdown.
%   [menu, suggested, allLibs] = sce_circ_enrichr_menu() fetches every Enrichr
%   library and returns MENU (cellstr) with a curated "Suggested" block on top
%   (newest GO BP/MF/CC, KEGG mouse, Reactome, WikiPathways mouse, MSigDB
%   Hallmark that actually exist), a separator, then all remaining libraries.
%
%   The GUI popupmenu uses MENU as its 'String'. Rows starting with the star or
%   the dashes are headers/separators — the callback should ignore them.
    d = webread('https://maayanlab.cloud/Enrichr/datasetStatistics');
    if ~isstruct(d); d = jsondecode(char(reshape(d, 1, []))); end
    st = d.statistics;
    if iscell(st)
        allLibs = string(cellfun(@(x) string(x.libraryName), st));
    else
        allLibs = string({st.libraryName});
    end
    allLibs = sort(allLibs(:));

    % preferred families, best-first; newest year of each is chosen
    prefs = ["GO_Biological_Process", "GO_Molecular_Function", "GO_Cellular_Component", ...
             "KEGG.*Mouse", "Reactome", "WikiPathway.*Mouse", "MSigDB_Hallmark"];
    suggested = strings(0, 1);
    for p = prefs
        hits = allLibs(~cellfun(@isempty, regexp(allLibs, p, 'once')));
        if ~isempty(hits); suggested(end+1, 1) = pick_newest(hits); end %#ok<AGROW>
    end
    suggested = unique(suggested, 'stable');
    rest      = setdiff(allLibs, suggested, 'stable');

    menu = [ {'★ ─ Suggested ─'}; cellstr(suggested(:)); ...
             {'──── all libraries ────'}; cellstr(rest(:)) ];
end

function s = pick_newest(hits)
% pick the entry with the largest 4-digit year (falls back to the last one)
    y = zeros(numel(hits), 1);
    for i = 1:numel(hits)
        tok = regexp(hits(i), '(\d{4})', 'tokens', 'once');
        if ~isempty(tok); y(i) = str2double(tok{1}); end
    end
    [~, idx] = max(y);
    s = hits(idx);
end
