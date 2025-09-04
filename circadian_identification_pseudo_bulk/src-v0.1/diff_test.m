function diff_bool = diff_test(fval,R0)
    
    % Perform difference of next point and current point of R0
    val0 = diff(R0);
    val0n = val0<0;
    val0p = val0>0;
    % Binarize
    val0 = val0p-val0n;
    
    % Perform difference of next point and current point of fval 
    val = diff(fval);
    val0n = val<0;
    val0p = val>0;
    val = val0p - val0n;

    valdiff = val - val0;

    % Test if change is same sign in end points
    nv = length(valdiff);
    idx = [1 nv];
    diff_bool0 = any(valdiff(idx));

    %diff_bool = any(valdiff);
    % Test if there are less than 30% of difference
    valdiff = sum(0.5*abs(valdiff))/nv;
    diff_bool1 = valdiff > 0.3;

    diff_bool = or(diff_bool0,diff_bool1) ;
end