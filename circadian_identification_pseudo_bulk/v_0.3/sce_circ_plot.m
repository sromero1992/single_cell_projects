function sce_circ_plot(sce, tmeta, cust_cells, plot_type, period12)
    % Plot identified circadian genes according to criterias
    % INPUT:
    % sce =========> Single cell experiment object
    % tmeta =======> Metatable containing circadian information (see usage)
    % period12 ====> Use period 12 (true) or 24 (false)
    % cust_cells ==> Custom cell type to compute
    % OUTPUT: 
    % Directories with plots in:
    % $(cust_cells)_circadian_syncronized_neighboars
    % $(cust_cells)_circadian_fit_expression_all
    % $(cust_cells)_circadian_fit_expression_failed
    % USAGE:
    % old_labels = ["1La" "2La" "3La" "4La" "5Da" "6Da" "7Da" "8Da"]';
    % new_labels = ["ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21"]';
    % times = [0: 3: 21]';
    % tmeta = table( old_labels, new_labels, times);
    % cust_genes = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
    %              "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
    % plot_type = 1;
    % sce_circ_plot(sce, tmeta, cust_cells, plot_type)

    if nargin < 4 || isempty(plot_type); plot_type = 1; end
    if nargin < 5 || isempty(period12); period12 = false; end

    % Test all genes always
    custom_genelist = [];
    % Required to input cell type or subcelltype
    rm_low_conf = true;

    % Define time variables
    t0 = tmeta.times(1);
    % This will work whenever times are evenly spaced
    tint = mean( diff(tmeta.times) );
    tf = tmeta.times(end);
    % This depends on different intervals and experimental points
    t = t0 : tint : tf;
    tval = t0 : 0.1 : tf;

    % Compute circadian information for cust_cells
    [T1, T2] = sce_circ_phase_estimation(sce, tmeta, rm_low_conf, period12, ...
                      custom_genelist, cust_cells );

    disp( "Final number of circadian genes: " + size(T1,1) )

    % Check if there are any classical circadian genes
    classic_circ = [ "arn" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
                     "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
                     "nfkb" "per" "wnt" "nrd" "rev"];
    
    ngenes = length(T1.Genes);
    ncirc = length(classic_circ);
    founds = false(ngenes,1);
    for ig = 1:ngenes
        for ic = 1:ncirc
            if startsWith( lower(T1.Genes(ig)), classic_circ(ic))
                founds(ig) = true;
                break;
            end
        end
    end
    
    if any(founds) 
        disp("Classic circadian genes identified:")
        disp(T1.Genes(founds))
    else
        disp("No classic circadian genes identified")
    end

    % Plot classic circadian genes and syncronized neighboars
    gidx = find(founds == 1)';
    colors = colormap(jet);
    path = strcat(cust_cells,'_circadian_syncronized_neighboars');
    disp("Syncronized gene results will be stored in "+path)
    mkdir(path);
    for i = gidx
        %disp(i)
        icolor = 1;
        ii = 1;
        f = figure('visible','off');
        gjdx = abs( T1.Acrophase_24(:) - T1.Acrophase_24(i) ) < 0.1;
        gjdx = find( gjdx == 1)';
        %disp(T2.Genes(gjdx));
        for j = gjdx
            Rzts = table2array( T2(j,2:end) );
            %disp(Rzts)
            plot(t, Rzts, 'Color', colors(icolor, :), 'Marker', 'o');
            icolor = icolor + 10;
            % Re-setting colors
            if icolor > 255; icolor = 1; end
            legend_text{ii} = T2.Genes(j);
            ii = ii + 1;
            hold on;
        end
        xlim([t0 tf])
        legend(legend_text, 'Location', 'bestoutside');
        title("Acrophase: "+T1.Acrophase_24(i));
        xlabel('Time (hrs)');
        ylabel('Expression');
        % Save figure
        fname = strcat("/plot_circadian_expr_sync_neig_", T1.Genes(i));
        fname = strcat(path,fname,".png");
        saveas(f,fname)
        clear legend_text;
        % Save text file with genes
        fname = strcat("/list_circadian_expr_sync_neig_", T1.Genes(i));
        fname = strcat(path,fname);
        tab = table(T2.Genes(gjdx));
        writetable(tab, fname,  'WriteVariableNames', 0);
        hold off;
    end

    switch plot_type
        case 1
            % Plot failed but confident circadian gene expression
            path = strcat(cust_cells,'_circadian_fit_expression_failed');
            mkdir(path);
            disp("Failed circadian gene plots will be stored in "+path)
            gidx = find(T1.Failed > 0)';
            for i = gidx
                f = figure('visible','off');
                fval = T1.Amp(i).*cos( 2*pi*( tval - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
                Rzts = table2array( T2(i,2:end) );
                plot(tval, fval);
                hold on;
                plot(t, Rzts);
                xlim([t0 tf])
                title("Gene - "+T1.Genes(i)+ " | pvalue " + ...
                      string( T1.pvalues(i) ) + " | failure: " + ... 
                      T1.Failed(i));
                xlabel('Time (hrs)');
                ylabel('Expression');
                legend({'Sine-fitted expression','Expression'},'Location','northwest');
                fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
                fname = strcat(path,fname,".png");
                saveas(f,fname)
                hold off;
            end
        case 2
            % Plot all confident circadian gene expression
            path = strcat(cust_cells,'_circadian_fit_expression_all');
            mkdir(path);
            disp("Results will be stored in "+path)
            for i = 1:ngenes
                f = figure('visible','off');
                fval = T1.Amp(i).*cos( 2*pi*( tval - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
                Rzts = table2array( T2(i,2:end) );
                plot(tval, fval);
                hold on;
                plot(t,Rzts);
                xlim([t0 tf])
                title("Gene - "+T1.Genes(i)+" | pvalue " + string(T1.pvalues(i)));
                xlabel('Time (hrs)');
                ylabel('Expression');
                legend({'Sine-fitted expression','Expression'},'Location','northwest');
                fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
                fname = strcat(path,fname,".png");
                saveas(f,fname)
                hold off;
            end
    end

end