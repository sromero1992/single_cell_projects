function [acrophase, amp, period, mesor, p_value] = estimate_phaseR_Ftest(Xg_zts, time_step, period12)
    % This function estimates the phase, amplitude, period, and mesor from gene expression data
    % using a sine function fit. It also computes the p-value for the significance of the sine fit 
    % using an F-test comparing it to a null model (mean model).
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
    %   - p_value: The p-value from the F-test indicating the significance of the sine model fit.

    % Number of time points
    nzts = size(Xg_zts, 2);

    % Initialize arrays to store the expression data and corresponding time points
    icells = cellfun(@length, Xg_zts); % Number of cells at each time point
    num_cells = sum(icells); % Total number of cells across all time points

    R = zeros(num_cells, 1); % Flattened array for expression data
    time_grid = zeros(num_cells, 1); % Corresponding time points for each cell

    % Reshape the data from cell array to flat arrays for fitting
    ic = 0;
    for it = 1:nzts
        % Assign expression data
        R(ic + 1:ic + icells(it)) = Xg_zts{it}(:);
        % Assign corresponding time points
        time_grid(ic + 1:ic + icells(it)) = (it - 1) * time_step;
        ic = ic + icells(it);
    end
    
    % Null model (mean model)
    mean_model = mean(R); % Mean of all expression data
    SSR_null = sum((R - mean_model).^2); % Sum of squared residuals for the null model

    % Define sine model based on the selected period (12 or 24 hours)
    if period12
        ft = fittype('amp * cos(2*pi*(t - acro)/12) + mesor', 'coefficients', {'acro', 'amp', 'mesor'}, 'independent', {'t'});
        period = 12;
    else
        ft = fittype('amp * cos(2*pi*(t - acro)/24) + mesor', 'coefficients', {'acro', 'amp', 'mesor'}, 'independent', {'t'});
        period = 24;
    end
    
    % Fit options for the sine model
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Algorithm', 'Trust-Region', ...
                         'Lower', [-24, -Inf, -Inf], ...
                         'Upper', [24, Inf, Inf], ...
                         'StartPoint', [0, mean_model * 2, mean_model]);
    
    % Fitting the sine model to the data
    [fmdl, gof] = fit(time_grid, R, ft, options);
    SSR_sine = gof.sse; % Sum of squared errors for the sine model

    % Calculate the F-statistic
    d1 = 2; % Difference in number of parameters (sine has 3, mean has 1)
    d2 = length(R) - 3; % Degrees of freedom for residuals (total observations minus sine model parameters)

    F_stat = ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2);

    % Calculate the p-value from the F-distribution
    p_value = 1 - fcdf(F_stat, d1, d2);
    
    % Outputs: the estimated parameters from the sine fit
    acrophase = fmdl.acro;
    amp = fmdl.amp;
    mesor = fmdl.mesor;
end