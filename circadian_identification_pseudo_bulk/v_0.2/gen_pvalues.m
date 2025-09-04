function [T1, confident_idx] = gen_pvalues( T1, mode, err_type)
    % This functions fits the errors from all time points to a distribution
    % log transform the errors to a nicer fit
    % INPUT:
    % T1       : contains the table having circadian gene predictions
    % mode = 1 : kernel function (Bayesian statistics)
    % mode = 2 : normal distribution
    % Keep_all : will keep all the results and will not filter pval < 0.05
    % OUTPUT:
    % T1 including pvalues column
    % keep_idx for keeping only confident fitted genes
    % USAGE: 
    % [T1, conf_idx] = gen_pvalues( T1, 1);
    % T1 = T1(keep_idx, :);

    switch err_type
        case "MAE_rel"
            func_val =  T1.MAE_rel; % MAE, RMSE, MAE_rel
        case "MAE"
            func_val =  T1.MAE;
    end

    log_func = log(func_val);

    switch mode
        case 1
            % Create a kernel distribution object by fitting it to the data.
            % Use the Epanechnikov kernel function (Bayesian statistics).
            pd_logMAE_skewnormal = fitdist(log_func, 'Kernel', 'Kernel', 'epanechnikov');

            % Find the threshold error value for p-value of 0.05 in the log-transformed space
            threshold_logMAE_skewnormal = icdf(pd_logMAE_skewnormal, 0.05);

            % Identify genes with MAE below the threshold in the original space
            significantGenesIdx_logMAE_skewnormal = func_val <= exp(threshold_logMAE_skewnormal);

            % compute cumulative distribution function to get pvalues for
            % all genes
            pvalues = cdf(pd_logMAE_skewnormal, log_func);

            % % Plot the histogram of log-transformed MAE and the fitted skew normal distribution
            % x_logMAE = linspace(min(log_func), max(log_func), 100);
            % y_logMAE = pdf(pd_logMAE_skewnormal, log_func);
            % figure;
            % histogram(log_func, 'Normalization', 'pdf');
            % hold on;
            % plot(x_logMAE, y_logMAE, 'r-', 'LineWidth', 2);
            % hold off;
            % title('Histogram of Log-transformed MAE and Fitted Skew Normal Distribution');
            % xlabel('Log-transformed Mean Absolute Error (Log(MAE))');
            % ylabel('Probability Density');

            % Remove non- significant pvalue predictions
            T1.pvalues = pvalues;
            confident_idx = significantGenesIdx_logMAE_skewnormal;

        case 2
            % % Fit a normal distribution to the log-transformed MAE data
            pd_logMAE_normal = fitdist(log_func, 'Normal');

            % Find the threshold error value for p-value of 0.05 in the log-transformed space
            threshold_logMAE_normal = icdf(pd_logMAE_normal, 0.05);

            % Identify genes with MAE below the threshold in the original space
            significantGenesIdx_logMAE_normal = func_val <= exp(threshold_logMAE_normal);

            % compute cumulative distribution function to get pvalues for
            % all genes
            pvalues = cdf(pd_logMAE_normal, log_func);

            % % Plot the histogram of log-transformed MAE and the fitted normal distribution
            % x_logMAE = linspace(min(log_MAE), max(log_MAE), 100);
            % y_logMAE = pdf(pd_logMAE_normal, x_logMAE);
            % figure;
            % histogram(log_MAE, 'Normalization', 'pdf');
            % hold on;
            % plot(x_logMAE, y_logMAE, 'r-', 'LineWidth', 2);
            % hold off;
            % title('Histogram of Log-transformed MAE and Fitted Normal Distribution');
            % xlabel('Log-transformed Mean Absolute Error (Log(MAE))');
            % ylabel('Probability Density');

            % Remove non- significant pvalue predictions
            T1.pvalues = pvalues;
            confident_idx = significantGenesIdx_logMAE_normal;
    end
end