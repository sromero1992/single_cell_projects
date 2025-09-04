function [T1, confident_idx] = gen_pvalues2(T1, mode)
    % This function fits the errors from all time points to a distribution
    % log transform the errors to a nicer fit
    % INPUT:
    % T1       : contains the table having circadian gene predictions
    % mode = 1 : kernel function (Bayesian statistics)
    % mode = 2 : normal distribution
    % OUTPUT:
    % T1 including combined pvalues column
    % confident_idx for keeping only confident fitted genes
    % USAGE: 
    % [T1, conf_idx] = gen_pvalues(T1, 1);
    % T1 = T1(keep_idx, :);

    % Define a nested function to calculate p-values
    function pvalues = calculate_pvalues(func_val, mode)
        log_func = log(func_val);
        switch mode
            case 1
                % Create a kernel distribution object by fitting it to the data.
                pd_logMAE_skewnormal = fitdist(log_func, 'Kernel', 'Kernel', 'epanechnikov');
                % Compute cumulative distribution function to get pvalues
                pvalues = cdf(pd_logMAE_skewnormal, log_func);
            case 2
                % Fit a normal distribution to the log-transformed data
                pd_logMAE_normal = fitdist(log_func, 'Normal');
                % Compute cumulative distribution function to get pvalues
                pvalues = cdf(pd_logMAE_normal, log_func);
        end
    end

    % Calculate p-values for MAE
    mae_pvalues = calculate_pvalues(T1.MAE, mode);
    T1.pvalues_mae = mae_pvalues;

    % Calculate p-values for MAE_rel
    mae_rel_pvalues = calculate_pvalues(T1.MAE_rel, mode);
    T1.pvalues_mae_rel = mae_rel_pvalues;

    % Combine p-values using Fisher's method
    combined_pvalues = -2 * (log(mae_pvalues) + log(mae_rel_pvalues)); % Chi2
    chi2_threshold = chi2inv(0.95, 4); % 95% confidence level, 4 degrees of freedom

    % Identify significant genes
    significantGenesIdx_combined = combined_pvalues <= chi2_threshold;

    % Add combined p-values to the table
    T1.pvalues_combined = chi2cdf(combined_pvalues, 4, 'upper'); % Upper tail to get p-values

    % Return confident gene indices
    confident_idx = significantGenesIdx_combined;
end