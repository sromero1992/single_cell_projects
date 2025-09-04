function sce_circ_plot(sce, tmeta, cust_cells, plot_type, period12, norm_str)
    % Plot identified circadian genes according to criteria
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
    if nargin < 6 || isempty(norm_str); norm_str ='lib_size' ; end

    % Define time variables
    tmeta.times = sortrows(tmeta.times);
    t0 = tmeta.times(1);
    % This should be improved for may times...
    tint = mean(diff(tmeta.times));
    
    disp("Time steps are : " + tint);
    tf = tmeta.times(end);
    t = t0 : tint : tf;
    tval = t0 : 0.1 : tf;
    
    cust_gene = [];
    % if plot_type == 3
    %     % Check if there are any classical circadian genes
    %     classic_circ = ["arnt" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
    %                     "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
    %                     "nfkb" "per" "wnt" "nr1d" "rev" "pik" "ror"];
    % 
    %     ngenes = length(sce.g);
    %     ncirc = length(classic_circ);
    %     idx = false(ngenes, 1);
    %     for ig = 1:ngenes
    %         for ic = 1:ncirc
    %             if startsWith( lower(sce.g(ig)), classic_circ(ic))
    %                 idx(ig) = true;
    %                 break;
    %             end
    %         end
    %     end
    %     cust_gene = sce.g(idx)';
    % end

    % Compute circadian information for cust_cells
    [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, plot_type == 1, ...
                                           period12, cust_gene, cust_cells, false, ...
                                           norm_str);

    disp("Final number of circadian genes: " + size(T1,1))

    % Check if there are any classical circadian genes
    classic_circ = ["arnt" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
                    "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
                    "nfkb" "per" "wnt" "nr1d" "rev" "pik" "ror"];

    ngenes = length(T1.Genes);
    ncirc = length(classic_circ);
    founds = false(ngenes,1);
    for ig = 1:ngenes
        for ic = 1:ncirc
            if startsWith(lower(T1.Genes(ig)), classic_circ(ic))
                founds(ig) = true;
                break;
            end
        end
    end

    % Define paths and file names
    if period12
        cust_cells = strcat(cust_cells, "_period_12_");
    else
        cust_cells = strcat(cust_cells, "_period_24_");
    end

    if any(founds) 
        disp("Classic circadian genes identified:")
        disp(T1.Genes(founds))
        fname = strcat(cust_cells, '_circadian_gene_list.txt');
        T1bak = T1(founds,:);
        writetable(T1bak, strcat(cust_cells, '_macro_circadian_analysis.csv'));
        writematrix(T1.Genes(founds), fname)
    else
        disp("No classic circadian genes identified")
    end

    % Extract necessary data for below cases
    gidx_conf = find(T1.pvalue <= 0.05)';
    gidx_nonconf = find(T1.pvalue > 0.05)';
    gidx_classic = find(founds == 1)';

    % Prepare paths
    path_sync = strcat(cust_cells, '_circadian_syncronized_phase');
    path_conf = strcat(cust_cells, '_circadian_confident');
    path_nonconf = strcat(cust_cells, '_circadian_non_confident');
    path_classic = strcat(cust_cells, '_classic_circadian');

    % Plot confident circadian genes
    if plot_type == 1
        mkdir(path_conf);
        disp("Confident circadian gene plots will be stored in " + path_conf);
        %parfor i = gidx_conf
        for i = gidx_conf
            f = figure('visible', 'off');
            amp = T1.Amp(i);
            acro = T1.Acrophase(i);
            T = T1.Period(i);
            mesor = T1.Mesor(i);
            fval = amp .* cos(2 * pi * (tval - acro) ./ T) + mesor;
            Rzts = table2array(T2(i, 2:end));
            plot(tval, fval);
            hold on;
            plot(t, Rzts);
            xlim([t0 tf])
            title("Gene - " + T1.Genes(i) + " | pvalue FT " + string(T1.pvalue(i)) + " | pvalue cor " + string(T1.pvalue_corr(i)));
            xlabel('Time (hrs)');
            ylabel('Expression');
            legend({'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
            fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
            fname = strcat(path_conf, fname, ".png");
            saveas(f, fname);
            hold off;
        end
    end

    % Plot non-confident circadian genes
    if plot_type == 2
        mkdir(path_nonconf);
        disp("Non-confident circadian gene plots will be stored in " + path_nonconf);
        %parfor i = gidx_nonconf
        for i = gidx_nonconf
            f = figure('visible', 'off');
            fval = T1.Amp(i) .* cos(2 * pi * (tval - T1.Acrophase(i)) ./ T1.Period(i)) + T1.Mesor(i);
            Rzts = table2array(T2(i, 2:end));
            plot(tval, fval);
            hold on;
            plot(t, Rzts);
            xlim([t0 tf])
            title("Gene - " + T1.Genes(i) + " | pvalue FT " + string(T1.pvalue(i)) + " | pvalue cor " + string(T1.pvalue_corr(i)));
            xlabel('Time (hrs)');
            ylabel('Expression');
            legend({'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
            fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
            fname = strcat(path_nonconf, fname, ".png");
            saveas(f, fname);
            hold off;
        end
    end

    % Plot classic circadian genes
    if plot_type == 3
        mkdir(path_classic);
        disp("Classical circadian gene results will be stored in " + path_classic);
        %parfor i = gidx_classic
        for i = gidx_classic
            f = figure('visible', 'off');
            fval = T1.Amp(i) .* cos(2 * pi * (tval - T1.Acrophase(i)) ./ T1.Period(i)) + T1.Mesor(i);
            Rzts = table2array(T2(i, 2:end));
            plot(tval, fval);
            hold on;
            plot(t, Rzts);
            xlim([t0 tf])
            title("Gene - " + T1.Genes(i) + " | pvalue FT " + string(T1.pvalue(i)) + " | pvalue cor " + string(T1.pvalue_corr(i)));
            xlabel('Time (hrs)');
            ylabel('Expression');
            legend({'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
            fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
            fname = strcat(path_classic, fname, ".png");
            saveas(f, fname);
            hold off;
        end
    end
    
    % Plot synchronized genes
    if plot_type == 4
        mkdir(path_sync);
        disp("Syncronized gene results will be stored in " + path_sync);
        for i = gidx_classic
            icolor = 1;
            ii = 1;
            f = figure('visible','off');
            gjdx = abs(T1.Acrophase_24 - T1.Acrophase_24(i)) < 0.01;
            gjdx = find(gjdx)';
            colors = colormap(jet);

            for j = gjdx
                Rzts = table2array(T2(j,2:end));
                plot(t, Rzts, 'Color', colors(icolor, :), 'Marker', 'o');
                icolor = icolor + 10;
                if icolor > 255; icolor = 1; end
                legend_text{ii} = T2.Genes(j);
                ii = ii + 1;
                hold on;
            end

            xlim([t0 tf])
            legend(legend_text, 'Location', 'bestoutside');
            title("Acrophase: " + T1.Acrophase_24(i));
            xlabel('Time (hrs)');
            ylabel('Expression');

            % Save figure
            fname = strcat("/plot_circadian_expr_sync_neig_", T1.Genes(i));
            fname = strcat(path_sync, fname, ".png");
            saveas(f, fname);
            clear legend_text;

            % Save text file with genes
            fname = strcat("/list_circadian_expr_sync_neig_", T1.Genes(i));
            fname = strcat(path_sync, fname);
            tab = table(T2.Genes(gjdx));
            writetable(tab, fname, 'WriteVariableNames', 0);
            hold off;
        end
    end
end