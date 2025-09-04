function sce_circ_plot_gene(sce, tmeta, cust_cells, plot_type, period12, cust_gene, handles)
    % Plot identified circadian genes according to criteria
    % INPUT:
    % sce =========> Single cell experiment object
    % tmeta =======> Metatable containing circadian information (see usage)
    % period12 ====> Use period 12 (true) or 24 (false)
    % cust_cells ==> Custom cell type to compute
    % handles =====> GUI handles to update the canvas with plots

    if nargin < 4 || isempty(plot_type); plot_type = 1; end
    if nargin < 5 || isempty(period12); period12 = false; end
    if nargin < 6 || isempty(cust_gene); cust_gene = []; end
    if nargin < 7; error('GUI handles are required.'); end

    % Required to input cell type or subcelltype
    rm_low_conf = false;

    % Define time variables
    t0 = tmeta.times(1);
    tint = mean(diff(tmeta.times));
    tf = tmeta.times(end);
    t = t0 : tint : tf;
    tval = t0 : 0.1 : tf;

    % Compute circadian information for cust_cells
    disp(tmeta)
    [T1, T2] = sce_circ_phase_estimation(sce, tmeta, rm_low_conf, period12, ...
                      cust_gene, cust_cells);

    % Plot classic circadian genes and synchronized neighbors
    gidx = find(founds == 1)';
    colors = colormap(jet);

    % Plot for each gene
    for i = gidx
        axes_handle = handles.PlotCanvas; % Adjust to your GUI axes handle
        cla(axes_handle); % Clear the previous plot
        
        icolor = 1;
        ii = 1;
        gjdx = abs(T1.Acrophase_24(:) - T1.Acrophase_24(i)) < 0.1;
        gjdx = find(gjdx == 1)';
        
        % Plot each synchronized neighbor
        for j = gjdx
            Rzts = table2array(T2(j,2:end));
            plot(axes_handle, t, Rzts, 'Color', colors(icolor, :), 'Marker', 'o');
            icolor = icolor + 10;
            if icolor > 255; icolor = 1; end
            legend_text{ii} = T2.Genes(j);
            ii = ii + 1;
            hold(axes_handle, 'on');
        end
        xlim(axes_handle, [t0 tf])
        legend(axes_handle, legend_text, 'Location', 'bestoutside');
        title(axes_handle, "Acrophase: "+T1.Acrophase_24(i));
        xlabel(axes_handle, 'Time (hrs)');
        ylabel(axes_handle, 'Expression');
        hold(axes_handle, 'off');
        
        % Update GUI canvas with the plot
        drawnow;
    end

    % Plot failed but confident circadian gene expression
    if plot_type == 1
        gidx = find(T1.Failed > 0)';
        for i = gidx
            axes_handle = handles.PlotCanvas; % Adjust to your GUI axes handle
            cla(axes_handle); % Clear the previous plot
            
            fval = T1.Amp(i).*cos(2*pi*(tval - T1.Acrophase(i))./T1.Period(i)) + T1.Mesor(i);
            Rzts = table2array(T2(i,2:end));
            plot(axes_handle, tval, fval);
            hold(axes_handle, 'on');
            plot(axes_handle, t, Rzts);
            xlim(axes_handle, [t0 tf])
            title(axes_handle, "Gene - "+T1.Genes(i)+ " | pvalue " + string(T1.pvalues(i)) + " | failure: " + T1.Failed(i));
            xlabel(axes_handle, 'Time (hrs)');
            ylabel(axes_handle, 'Expression');
            legend(axes_handle, {'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
            hold(axes_handle, 'off');
            
            % Update GUI canvas with the plot
            drawnow;
        end
    end

    % Plot all confident circadian gene expression
    if plot_type == 2
        for i = 1:length(T1.Genes)
            axes_handle = handles.PlotCanvas; % Adjust to your GUI axes handle
            cla(axes_handle); % Clear the previous plot
            
            fval = T1.Amp(i).*cos(2*pi*(tval - T1.Acrophase(i))./T1.Period(i)) + T1.Mesor(i);
            Rzts = table2array(T2(i,2:end));
            plot(axes_handle, tval, fval);
            hold(axes_handle, 'on');
            plot(axes_handle, t, Rzts);
            xlim(axes_handle, [t0 tf])
            title(axes_handle, "Gene - "+T1.Genes(i)+" | pvalue " + string(T1.pvalues(i)));
            xlabel(axes_handle, 'Time (hrs)');
            ylabel(axes_handle, 'Expression');
            legend(axes_handle, {'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
            hold(axes_handle, 'off');
            
            % Update GUI canvas with the plot
            drawnow;
        end
    end
end