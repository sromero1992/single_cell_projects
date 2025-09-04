function [acrophase, amp, period, mesor] = estimate_phaseR(R, ...
                                                    time_cycle, time_step, period12)
    % sce_in contains the desired cell type SCE
    % gs is the list of genes to look at on SCE
    % time is the time frame within batches 

    % Time grid for fitting 
    n_bulk = size(R,1);
    nzts = size(R,2);
    batch_grid = 1:sum(nzts);
    time_grid = (batch_grid-1)*time_step;
    time_grid = time_grid';
    time_grid = repmat(time_grid, n_bulk, 1);
    medr = mean(R, "all");
    % Re-shaping R from 2D to 1D 
    R = reshape(R',[],1);

    % Sine function fit
    if period12
        ft = fittype('amp * cos( 2*pi*(t - acro)/12) + mesor','coefficients',...
                    {'acro','amp','mesor'}, 'independent', {'t'});
        options = fitoptions('Method','NonlinearLeastSquares',...
                  'Algorithm', 'Trust-Region',...
                  'Lower', [-24 -Inf -Inf ], ...
                  'Upper', [24 Inf Inf],...
                  'StartPoint',[0, medr*2, medr]);
        % Fitting function here
        [fmdl, ~] = fit(time_grid, R, ft, options);
        acrophase = fmdl.acro;
        amp = fmdl.amp;
        mesor = fmdl.mesor;
        period = 12;
    else
        ft = fittype('amp * cos( 2*pi*(t - acro)/24) + mesor','coefficients',...
                    {'acro','amp','mesor'}, 'independent', {'t'});
        options = fitoptions('Method','NonlinearLeastSquares',...
                  'Algorithm', 'Trust-Region',...
                  'Lower', [-24 -Inf -Inf ], ...
                  'Upper', [24 Inf Inf],...
                  'StartPoint',[0, medr*2, medr]);
        % Fitting function here
        [fmdl, ~] = fit(time_grid, R, ft, options);
        acrophase = fmdl.acro;
        amp = fmdl.amp;
        mesor = fmdl.mesor;
        period = 24;
    end

end