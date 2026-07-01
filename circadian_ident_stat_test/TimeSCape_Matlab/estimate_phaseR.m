function [acrophase, amp, period, mesor, ...
           p_value, rho, p_value_macro] = estimate_phaseR(Xg_zts, actual_times, ...
            period12, test_type)
% ESTIMATE_PHASER
%   Fits a cosine (cosinor) model to single-cell gene expression data and
%   tests its significance against an intercept-only null model.
%
%   The function accepts ACTUAL ZT times for each slot, so it works
%   correctly even when one or more time points are missing from the data
%   (the cosine is fitted at the true hour values, not assumed-uniform
%   positions).
%
% INPUTS:
%   Xg_zts      - Cell array (1 × nzts). Each cell holds the expression
%                 values of every cell collected at that time point.
%                 Empty cells (no cells at that time) are skipped
%                 automatically.
%   actual_times - Numeric vector (1 × nzts) of the real ZT hours
%                  corresponding to each slot in Xg_zts, e.g.
%                  [0, 3, 6, 9, 15, 18, 21] for 7 evenly-spaced points
%                  with ZT12 absent.
%   period12    - Logical; true → 12-hr period, false → 24-hr period.
%   test_type   - 'Ftest' or 'LRT'.
%
% OUTPUTS:
%   acrophase    - Estimated time of peak expression (hrs).
%   amp          - Estimated cosine amplitude.
%   period       - Period used (12 or 24 hrs).
%   mesor        - Midline estimating statistic of rhythm.
%   p_value      - F-test (or LRT) p-value.
%   rho          - Pearson correlation of cosine fit vs. R0 (mean/timepoint).
%   p_value_macro- p-value for that Pearson correlation.

    nzts = numel(Xg_zts);

    % ── Remove empty time slots ───────────────────────────────────────────
    % Keeps the function robust when a time point has no cells at all.
    non_empty = cellfun(@(x) ~isempty(x) && numel(x) > 0, Xg_zts);
    Xg_zts      = Xg_zts(non_empty);
    actual_times = actual_times(non_empty);
    nzts         = sum(non_empty);

    % Need at least 4 time points to fit 3-parameter cosine + meaningful df
    if nzts < 4
        acrophase    = NaN; amp   = NaN;
        period       = NaN; mesor = NaN;
        p_value      = NaN; rho   = NaN;
        p_value_macro = NaN;
        return;
    end

    % ── Flatten cells → single arrays ────────────────────────────────────
    icells    = cellfun(@numel, Xg_zts);
    num_cells = sum(icells);
    R         = zeros(num_cells, 1);
    time_grid = zeros(num_cells, 1);
    R0        = zeros(1, nzts);

    ic          = 0;
    max_R0      = -Inf;
    max_peak_t  = actual_times(1);

    for it = 1:nzts
        n_it = icells(it);
        R(ic+1 : ic+n_it)         = Xg_zts{it}(:);
        meanval                    = mean(Xg_zts{it}(:), 'omitnan');
        R0(it)                     = meanval;
        time_grid(ic+1 : ic+n_it) = actual_times(it);   % ← actual hour
        ic = ic + n_it;
        if meanval > max_R0
            max_R0     = meanval;
            max_peak_t = actual_times(it);
        end
    end

    % ── Null model ───────────────────────────────────────────────────────
    mean_model = mean(R);
    SSR_null   = sum((R - mean_model).^2);

    % ── Cosine model ─────────────────────────────────────────────────────
    if period12; period = 12; else; period = 24; end

    max_amp_guess = max_R0 - mean_model;

    ft = fittype(sprintf('amp * cos(2*pi*(t - acro)/%d) + mesor', period), ...
                 'coefficients', {'acro', 'amp', 'mesor'}, ...
                 'independent',  {'t'});

    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Algorithm',  'Trust-Region', ...
                         'Lower',      [0,      -Inf, 0  ], ...
                         'Upper',      [period,  Inf, Inf], ...
                         'StartPoint', [max_peak_t, max_amp_guess, mean_model]);

    [fmdl, gof] = fit(time_grid, R, ft, options);
    SSR_sine    = gof.sse;

    % ── Statistical test ─────────────────────────────────────────────────
    if strcmp(test_type, 'Ftest')
        d1      = 2;               % extra params in cosine vs. null
        d2      = num_cells - 3;   % residual df
        F_stat  = ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2);
        p_value = 1 - fcdf(F_stat, d1, d2);

    elseif strcmp(test_type, 'LRT')
        s2_null = SSR_null / num_cells;
        logL_null = -0.5 * num_cells * (log(2*pi*s2_null) + 1);
        s2_sine = SSR_sine / num_cells;
        logL_sine = -0.5 * num_cells * (log(2*pi*s2_sine) + 1);
        p_value = 1 - chi2cdf(-2*(logL_null - logL_sine), 2);

    else
        error('test_type must be "Ftest" or "LRT".');
    end

    acrophase = fmdl.acro;
    amp       = fmdl.amp;
    mesor     = fmdl.mesor;

    % ── Canonical form: positive amplitude ───────────────────────────────
    % If the NLS converges with amp < 0 the trough is at acrophase and the
    % true peak is half a period away. Flip both so the reported amplitude
    % is always positive and acrophase is the true peak time.
    if amp < 0
        amp       = -amp;
        acrophase = mod(acrophase + period/2, period);
    end
    acrophase = mod(acrophase, period);   % robust wrap into [0, period)

    % ── Correlation t-test: cosine fit vs. per-timepoint means ───────────
    % Uses actual ZT values so missing points don't bias the correlation.
    fval = amp * cos(2*pi*(actual_times - acrophase) / period) + mesor;
    [rho, p_value_macro] = corr(R0(:), fval(:), 'Type', 'Pearson');

end
