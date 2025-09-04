function [acrophase, amp, period, mesor, ... 
           p_value, rho, p_value_macro] = estimate_phaseR(Xg_zts, time_step, ...
            period12, test_type)
    % This function estimates the phase, amplitude, period, and mesor from gene expression data
    % using a sine function fit. It computes the p-value for the significance of the sine fit 
    % using either a Likelihood Ratio Test (LRT) or an F-test, as specified by the test_type argument.
    %
    % Inputs:
    %   - Xg_zts: A cell array where each cell contains the expression data of all cells at each time point.
    %   - time_step: The time difference between successive time points.
    %   - R0:      Mean values per time point from all cells.
    %   - period12: A boolean flag indicating whether to fit a 12-hour (true) or 24-hour (false) period sine function.
    %   - test_type: A string specifying the test type, either 'LRT' for Likelihood Ratio Test or 'Ftest' for F-test.
    %
    % Outputs:
    %   - acrophase: The estimated acrophase (time of peak expression).
    %   - amp: The estimated amplitude of the sine function.
    %   - period: The period used for fitting (12 or 24 hours).
    %   - mesor: The estimated mesor (mean level of expression).
    %   - p_value: The p-value from the specified test indicating the significance of the sine model fit.

    % Number of time points
    nzts = size(Xg_zts, 2);

    % Initialize arrays to store the expression data and corresponding time points
    icells = cellfun(@length, Xg_zts); % Number of cells at each time point
    num_cells = sum(icells); % Total number of cells across all time points
    R = zeros(num_cells, 1); % Flattened array for expression data
    time_grid = zeros(num_cells, 1); % Corresponding time points for each cell
    R0 = zeros(1, nzts); % Mean values per time point

    % Reshape the data from cell array to flat arrays for fitting
    ic = 0;
    max_amp = -1000;
    max_amp_old = -1000;
    max_peak_t = 0;
    for it = 1:nzts
        % Assign expression data
        R(ic + 1:ic + icells(it)) = Xg_zts{it}(:);
        meanval = mean( Xg_zts{it}(:), 'omitnan');
        R0(it) = meanval;
        max_amp = max( max_amp, meanval);
        % Assign corresponding time points
        time_grid(ic + 1:ic + icells(it)) = (it - 1) * time_step;
        ic = ic + icells(it);
        if max_amp > max_amp_old
            max_amp_old = max_amp;
            max_peak_t = (it - 1) * time_step;
        end
    end
    
    % Null model (mean model)
    mean_model = mean(R); % Mean of all expression data
    SSR_null = sum((R - mean_model).^2); % Sum of squared residuals for the null model

    % Guess amplitude
    max_amp_guess = max_amp - mean_model;

    %fprintf("Max peak %f mean %f and amp %f \n", max_peak_t, mean_model, max_amp)

    % Define sine model based on the selected period (12 or 24 hours)
    if period12
        ft = fittype('amp * cos(2*pi*(t - acro)/12) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 12;
        % Fit options for the sine model
        options = fitoptions('Method', 'NonlinearLeastSquares', ...
                             'Algorithm', 'Trust-Region', ...
                             'Lower', [0, -Inf, 0], ...
                             'Upper', [12, Inf, Inf], ...
                             'StartPoint', [max_peak_t, max_amp_guess, mean_model]);
    else
        ft = fittype('amp * cos(2*pi*(t - acro)/24) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 24;
        % Fit options for the sine model
        options = fitoptions('Method', 'NonlinearLeastSquares', ...
                             'Algorithm', 'Trust-Region', ...
                             'Lower', [0, -Inf, 0], ...
                             'Upper', [24, Inf, Inf], ...
                             'StartPoint', [max_peak_t, max_amp_guess, mean_model]);
    end
    
    % Fitting the sine model to the data
    [fmdl, gof] = fit(time_grid, R, ft, options);
    SSR_sine = gof.sse; % Sum of squared errors for the sine model

    if strcmp(test_type, 'Ftest')
        % F-test
        d1 = 2; % Difference in number of parameters (sine has 3, mean has 1)
        d2 = length(R) - 3; % Degrees of freedom for residuals (total observations minus sine model parameters)
        F_stat = ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2);
        p_value = 1 - fcdf(F_stat, d1, d2);
        
    elseif strcmp(test_type, 'LRT')
        % Likelihood Ratio Test (LRT)
        sigma2_null = SSR_null / num_cells;
        logL_null = -0.5 * num_cells * (log(2 * pi * sigma2_null) + 1);
        
        sigma2_sine = SSR_sine / num_cells;
        logL_sine = -0.5 * num_cells * (log(2 * pi * sigma2_sine) + 1);
        
        LRT_stat = -2 * (logL_null - logL_sine);
        df = 2; % Degrees of freedom difference between the sine and null model
        p_value = 1 - chi2cdf(LRT_stat, df);
    else
        error('Invalid test type. Choose either "LRT" or "Ftest".');
    end
   
    % Outputs: the estimated parameters from the sine fit
    acrophase = fmdl.acro;
    amp = fmdl.amp;
    mesor = fmdl.mesor;
    
    % Mean expression p-value test through correlation t-test
    t_pts = (0:nzts-1) * time_step;
    fval = amp * cos(2*pi*(t_pts - acrophase)/24) + mesor;

    [rho, p_value_macro] = corr(R0', fval', 'Type', 'Pearson'); 

end