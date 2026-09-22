function bins = sce_circ_bin_genes(T1_conf, bin_width, period)
%SCE_CIRC_BIN_GENES  Bin confident circadian genes by acrophase into ZT windows.
%   bins = sce_circ_bin_genes(T1_conf, bin_width) groups the genes in T1_conf by
%   their Acrophase_24 into windows of BIN_WIDTH hours (default 3 = one ZT
%   interval for an 8-point design). The ZT window is recorded per bin.
%
%   T1_conf   - table with columns Genes and Acrophase_24
%               (from *_circadian_analysis_confident.csv, read with readtable).
%   bin_width - window width in hours (default 3, user-definable).
%   period    - cycle length (default 24).
%
%   Returns a struct array with fields: label ('ZT06-09'), lo, hi, genes (cellstr).
    if nargin < 2 || isempty(bin_width); bin_width = 3;  end
    if nargin < 3 || isempty(period);    period    = 24; end
    edges = 0:bin_width:period;
    nb    = numel(edges) - 1;
    acro  = T1_conf.Acrophase_24;
    genes = string(T1_conf.Genes);
    bins  = struct('label', {}, 'lo', {}, 'hi', {}, 'genes', {});
    for i = 1:nb
        lo = edges(i); hi = edges(i+1);
        if i < nb
            sel = acro >= lo & acro <  hi;
        else
            sel = acro >= lo & acro <= hi;   % last bin closed on the right
        end
        bins(i).label = sprintf('ZT%02.0f-%02.0f', lo, hi);
        bins(i).lo    = lo;
        bins(i).hi    = hi;
        bins(i).genes = cellstr(genes(sel));
    end
end
