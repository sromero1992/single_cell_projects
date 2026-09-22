function libs = sce_circ_enrichr_libs(pattern)
%SCE_CIRC_ENRICHR_LIBS  List available Enrichr gene-set libraries (exact names).
%   libs = sce_circ_enrichr_libs()          returns all library names.
%   libs = sce_circ_enrichr_libs('GO')      returns names containing 'GO'.
%   libs = sce_circ_enrichr_libs('Mouse')   returns names containing 'Mouse'.
%   Use this to get the EXACT string to pass as the library to sce_circ_enrichr.
    d = webread('https://maayanlab.cloud/Enrichr/datasetStatistics');
    if ~isstruct(d); d = jsondecode(char(reshape(d, 1, []))); end
    st = d.statistics;
    if iscell(st)
        libs = string(cellfun(@(x) string(x.libraryName), st));
    else
        libs = string({st.libraryName});
    end
    if nargin >= 1 && ~isempty(pattern)
        libs = libs(contains(libs, pattern, 'IgnoreCase', true));
    end
    libs = sort(libs(:));
    disp(libs);
end
