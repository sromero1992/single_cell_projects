function T = sce_circ_enrichr(genes, library, adjp_cutoff, n_top)
%SCE_CIRC_ENRICHR  Over-representation analysis via the Enrichr web API.
%   T = sce_circ_enrichr(genes, library) submits GENES to Enrichr and returns a
%   table of enriched terms from LIBRARY, sorted by p-value. The per-term
%   "Genes" column holds the OVERLAP genes (input genes that are in the term),
%   i.e. the phase-restricted set to score with AUCell.
%
%   genes       - string/cellstr of gene symbols (one phase bin).
%   library     - Enrichr library, e.g. 'KEGG_2019_Mouse', 'Reactome_2022',
%                 'GO_Biological_Process_2021' (default 'KEGG_2019_Mouse').
%   adjp_cutoff - keep terms with adjusted p <= this (default 0.05).
%   n_top       - keep at most this many terms (default 10).
%
%   Requires internet access (Enrichr REST API, maayanlab.cloud).
    if nargin < 2 || isempty(library);     library     = 'KEGG_2019_Mouse'; end
    if nargin < 3 || isempty(adjp_cutoff); adjp_cutoff = 0.05; end
    if nargin < 4 || isempty(n_top);       n_top       = 10;   end
    import matlab.net.http.*
    import matlab.net.http.io.*

    genes = cellstr(string(genes));
    genes = genes(~cellfun(@isempty, genes));
    empty = cell2table(cell(0,5), 'VariableNames', {'Term','Pvalue','AdjPvalue','nOverlap','Genes'});
    if isempty(genes); T = empty; return; end
    listStr = strjoin(genes, newline);
    base    = 'https://maayanlab.cloud/Enrichr';

    % ---- addList (multipart/form-data POST) ----
    uri  = matlab.net.URI([base '/addList']);
    resp = [];
    for attempt = 1:6                                   % retry on 429 (rate limit)
        provider = MultipartFormProvider('list', listStr, 'description', 'TimeSCape');
        req = RequestMessage('POST', [], provider);
        try
            resp = req.send(uri);
        catch
            pause(1.5 * attempt); continue;
        end
        if resp.StatusCode == matlab.net.http.StatusCode.OK; break; end
        if double(resp.StatusCode) == 429; pause(1.5 * attempt); continue; end
        break;
    end
    if isempty(resp) || resp.StatusCode ~= matlab.net.http.StatusCode.OK
        sc = ''; if ~isempty(resp); sc = char(string(resp.StatusCode)); end
        error('Enrichr addList failed (%s) after retries. Wait a moment and try again.', sc);
    end
    raw = resp.Body.Data;                      % may be a raw JSON string
    if ischar(raw) || isstring(raw)
        payload = jsondecode(char(raw));
    else
        payload = raw;
    end
    userListId = payload.userListId;

    % ---- enrich (GET) ----
    url  = sprintf('%s/enrich?userListId=%d&backgroundType=%s', base, double(userListId), library);
    data = webread(url);
    if ~isstruct(data)
        txt = char(reshape(data, 1, []));      % force a ROW char vector for jsondecode
        try
            data = jsondecode(txt);
        catch
            error(['Enrichr returned non-JSON for library "%s" (name may be invalid). ' ...
                   'Run  sce_circ_enrichr_libs(''GO'')  to find the exact library name.'], library);
        end
    end
    fn   = fieldnames(data);
    if isempty(fn) || isempty(data.(fn{1})); T = empty; return; end
    rows = data.(fn{1});                 % cell array; each element is one term row
    n    = numel(rows);
    Term = strings(n,1); Pv = zeros(n,1); Adj = zeros(n,1); nOv = zeros(n,1); Gs = cell(n,1);
    for i = 1:n
        r = rows{i};
        Term(i) = string(r{2});
        Pv(i)   = r{3};
        Adj(i)  = r{7};
        Gs{i}   = cellstr(string(r{6}));   % overlap genes
        nOv(i)  = numel(Gs{i});
    end
    T = table(Term, Pv, Adj, nOv, Gs, 'VariableNames', {'Term','Pvalue','AdjPvalue','nOverlap','Genes'});
    T = T(T.AdjPvalue <= adjp_cutoff, :);
    T = sortrows(T, 'Pvalue', 'ascend');
    if height(T) > n_top; T = T(1:n_top, :); end
end
