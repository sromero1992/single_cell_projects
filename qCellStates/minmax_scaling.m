function X = minmax_scaling(X)
    % Scales X matrix to an interval from 0 to 1
   
    maxval = max(X,[],'all');
    minval = min(X,[],'all');

    X = (X - minval)./(maxval - minval);

end