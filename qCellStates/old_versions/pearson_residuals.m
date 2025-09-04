function z = pearson_residuals(X, theta)
    %{
    This function computes Pearson residual (put equation)
    INPUT:
    X -----> count matrix (g,c)
    theta -> parameter (100)
    OUTPU:
    z -----> Pearson residual matrix (g,c)
    %}
    if nargin < 2
        theta = 100;
    end
    X(isnan(X)) = 0;

    mu = (sum(X, 2) * sum(X, 1)) ./ sum(X(:));
    
    % NB variance
    sigma = sqrt(mu + mu.^2 / theta);

    z = (X - mu) ./ sigma;
    z(isnan(z)) = 0;

    n = size(X, 2);
    sn = sqrt(n);
    z(z >  sn) =  sn;
    z(z < -sn) = -sn;

end