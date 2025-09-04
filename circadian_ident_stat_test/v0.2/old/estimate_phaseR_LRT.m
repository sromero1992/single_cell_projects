function [acrophase, amp, period, mesor, p_value] = estimate_phaseR_LRT(Xg_zts, time_step, period12)
    % This function estimates the acrophase, amplitude, period, and mesor from gene expression data
    % using a sine function fit. It also computes the p-value for the significance of the sine fit
    % using a Likelihood Ratio Test (LRT) comparing it to a null model (mean model).
    %
    % Inputs:
    %   - Xg_zts: A cell array where each cell contains the expression data of all cells at each time point.
    %   - time_step: The time difference between successive time points.
    %   - period12: A boolean flag indicating whether to fit a 12-hour (true) or 24-hour (false) period sine function.
    %
    % Outputs:
    %   - acrophase: The estimated acrophase (time of peak expression).
    %   - amp: The estimated amplitude of the sine function.
    %   - period: The period used for fitting (12 or 24 hours).
    %   - mesor: The estimated mesor (mean level of expression).
    %   - p_value: The p-value from the Likelihood Ratio Test (LRT) indicating the significance of the sine model fit.

    % Number of time points (nzts) in the data
    nzts = size(Xg_zts, 2);

    % Get the number of cells at each time point
    icells = cellfun(@length, Xg_zts, 'UniformOutput', true);

    % Total number of cells across all time points
    num_cells = sum(icells);

    % Initialize arrays to hold expression data (R) and corresponding time points (time_grid)
    R = zeros(num_cells, 1);
    time_grid = zeros(num_cells, 1);

    % Flatten the data from cell array to vectors for fitting
    ic = 0;
    for it = 1:nzts
        % Assign expression data from the current time point
        R(ic + 1:ic + icells(it)) = Xg_zts{it}(:);
        % Assign corresponding time points
        time_grid(ic + 1:ic + icells(it)) = (it - 1) * time_step;
        % Update the index for the next batch of cells
        ic = ic + icells(it);
    end
    
    % Null model (mean model) - calculates the mean of all expression data
    mean_model = mean(R);
    
    % Sum of squared residuals (SSR) for the null model
    SSR_null = sum((R - mean_model).^2);

    % Estimate the variance of the noise under the null model
    sigma2_null = SSR_null / num_cells;

    % Corrected log-likelihood for the null model
    logL_null = -0.5 * num_cells * (log(2 * pi * sigma2_null) + 1);

    % Define the sine model with a period of either 12 or 24 hours
    if period12
        % 12-hour period sine model
        ft = fittype('amp * cos(2*pi*(t - acro)/12) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 12;
    else
        % 24-hour period sine model
        ft = fittype('amp * cos(2*pi*(t - acro)/24) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 24;
    end
    
    % Fit options for the sine model
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Algorithm', 'Trust-Region', ...
                         'Lower', [-24, -Inf, -Inf], ...
                         'Upper', [24, Inf, Inf], ...
                         'StartPoint', [0, mean_model * 2, mean_model]);
    
    % Fit the sine model to the data
    [fmdl, gof] = fit(time_grid, R, ft, options);

    % Sum of squared residuals (SSR) for the sine model
    SSR_sine = gof.sse;

    % Estimate the variance of the noise under the sine model
    sigma2_sine = SSR_sine / num_cells;

    % Corrected log-likelihood for the sine model
    logL_sine = -0.5 * num_cells * (log(2 * pi * sigma2_sine) + 1);

    % Calculate the Likelihood Ratio Test (LRT) statistic
    LRT_stat = -2 * (logL_null - logL_sine);

    % Degrees of freedom difference between the sine model (3 parameters) and the null model (1 parameter)
    df = 3 - 1;

    % Compute the p-value using the chi-square distribution
    p_value = 1 - chi2cdf(LRT_stat, df);
    
    % Output the estimated parameters from the sine fit
    acrophase = fmdl.acro;
    amp = fmdl.amp;
    mesor = fmdl.mesor;
end