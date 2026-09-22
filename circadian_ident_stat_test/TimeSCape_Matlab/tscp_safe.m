function s = tscp_safe(x)
%TSCP_SAFE  Canonical cell-type -> filesystem name (matches the analysis writer).
    s = regexprep(strtrim(char(x)), '[^a-zA-Z0-9_]', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
end
