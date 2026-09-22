function stats = sce_circ_pathway_plot(sce, tmeta, pathway_genes, pathway_name, cust_cells, period12, axHandle, color_name, use_violin, standardize, norm_str)
%SCE_CIRC_PATHWAY_PLOT  Circadian test + plot for a pathway's AUCell activity.
%   Scores cells for PATHWAY_GENES (AUCell), groups by ZT, fits a cosine via
%   estimate_phaseR, and plots per-ZT distributions + means + cosine on a white
%   grid matching the gene plot. Returns STATS (acrophase, amp, mesor, period,
%   pvalue, corr, pvalue_corr).
%
%   color_name  - single violin colour ('light blue','teal','orange',...).
%   use_violin  - true = violin, false = box.
%   standardize - true = z-score the AUCell score across cells (centred at 0);
%                 false = raw AUCell activity (>= 0).
%   norm_str    - 'lib_size' (default) or 'none'.
    if nargin < 6  || isempty(period12);    period12    = false;        end
    if nargin < 7  || isempty(axHandle);    figure; axHandle = axes;    end
    if nargin < 8  || isempty(color_name);  color_name  = 'light blue'; end
    if nargin < 9  || isempty(use_violin);  use_violin  = true;         end
    if nargin < 10 || isempty(standardize); standardize = false;        end
    if nargin < 11 || isempty(norm_str);    norm_str    = 'lib_size';   end
    if period12; period = 12; else; period = 24; end

    ic0 = find(sce.c_cell_type_tx == cust_cells);
    if isempty(ic0); error('No cells of type "%s".', cust_cells); end
    Xsub      = sce.X(:, ic0);
    batch_sub = sce.c_batch_id(ic0);
    if strcmp(norm_str, 'lib_size'); Xn = log1p(pkg.norm_libsize(Xsub, 1e4));
    else;                            Xn = Xsub; end

    score = sce_circ_aucell(Xn, sce.g, pathway_genes);
    if standardize
        mu = mean(score,'omitnan'); sd = std(score,'omitnan');
        if sd == 0 || isnan(sd); sd = 1; end
        score = (score - mu) / sd;
        ylab = 'Pathway activity (z-score)';
    else
        ylab = 'AUCell pathway activity';
    end

    batch_time = unique(batch_sub); nz = numel(batch_time); at = nan(nz,1);
    for it = 1:nz
        k = find(tmeta.new_labels == batch_time(it), 1);
        if ~isempty(k); at(it) = tmeta.ZT_times(k); end
    end
    valid = ~isnan(at); batch_time = batch_time(valid); at = at(valid);
    [at, ord] = sort(at); batch_time = batch_time(ord); nz = numel(at);

    Xg = cell(1, nz); means = zeros(1, nz);
    for it = 1:nz
        ics = find(batch_sub == batch_time(it));
        Xg{it}    = score(ics);
        means(it) = mean(score(ics), 'omitnan');
    end

    [acro, amp, per, mesor, pval, rho, pval_corr] = estimate_phaseR(Xg, at', period12, 'Ftest');
    stats = struct('acrophase',acro,'amp',amp,'mesor',mesor,'period',per, ...
                   'pvalue',pval,'corr',rho,'pvalue_corr',pval_corr);

    % ---- white theme + grid (match the gene plot) ----
    cla(axHandle, 'reset'); hold(axHandle, 'on');
    set(axHandle, 'Color',[1 1 1], 'XColor',[0.15 0.15 0.15], 'YColor',[0.15 0.15 0.15], ...
        'GridColor',[0.82 0.82 0.82], 'GridAlpha',0.7, 'Box','on');
    grid(axHandle, 'on');

    col   = local_color(color_name);
    halfw = period / (6 * max(nz,1));
    for it = 1:nz
        v = Xg{it}(:); v = v(~isnan(v));
        if use_violin && numel(v) > 5 && (max(v) > min(v))
            xi = linspace(min(v), max(v), 100);          % clip density to data range
            f  = ksdensity(v, xi);
            f  = f / max(f) * halfw;
            fill(axHandle, at(it) + [f, -fliplr(f)], [xi, fliplr(xi)], col, ...
                 'FaceAlpha', 0.55, 'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);
        elseif ~isempty(v)
            q = quantile(v, [0.25 0.5 0.75]);
            fill(axHandle, at(it)+[-halfw -halfw halfw halfw], [q(1) q(3) q(3) q(1)], ...
                 col, 'FaceAlpha',0.55, 'EdgeColor',[0.15 0.15 0.15], 'LineWidth',0.5);
            plot(axHandle, at(it)+[-halfw halfw], [q(2) q(2)], 'k-', 'LineWidth', 1);
        end
    end
    tt = min(at):0.1:max(at);
    ff = amp * cos(2*pi*(tt - acro)/period) + mesor;
    plot(axHandle, tt, ff, '-', 'Color',[0.13 0.47 0.71], 'LineWidth', 3);
    plot(axHandle, at, means, 'o-', 'Color',[0.85 0.10 0.10], 'LineWidth', 2, ...
         'MarkerFaceColor',[0.85 0.10 0.10]);
    acx = mod(acro, period);
    xline(axHandle, acx, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',1);
    yl = get(axHandle, 'YLim');
    text(axHandle, acx, yl(1) + 0.60*(yl(2)-yl(1)), sprintf('acro %.1f h', acx), ...
         'Color',[0.30 0.30 0.30], 'Rotation',90, 'FontSize',9, 'FontWeight','bold', ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
    hold(axHandle, 'off');
    xlabel(axHandle, 'ZT (h)', 'Color',[0.1 0.1 0.1]);
    ylabel(axHandle, ylab,     'Color',[0.1 0.1 0.1]);
    ttl = char(pathway_name);
    if numel(ttl) > 58; ttl = [ttl(1:55) '...']; end
    title(axHandle, sprintf('%s  (F p=%.2g, corr p=%.2g)', ttl, pval, pval_corr), ...
          'Interpreter','none', 'Color',[0.1 0.1 0.1], 'FontSize',10);
    set(axHandle, 'XTick', at);
end

function c = local_color(name)
    switch lower(char(name))
        case {'light blue','lightblue','blue'}; c = [0.55 0.75 0.95];
        case 'teal';   c = [0.20 0.60 0.60];
        case 'orange'; c = [0.95 0.60 0.25];
        case 'purple'; c = [0.55 0.40 0.75];
        case 'green';  c = [0.35 0.70 0.45];
        case 'gray';   c = [0.65 0.65 0.65];
        case 'red';    c = [0.85 0.45 0.45];
        otherwise;     c = [0.55 0.75 0.95];
    end
end
