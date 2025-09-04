%old_labels = ["1La" "2La" "3La" "4La" "5Da" "6Da" "7Da" "8Da"  ]';
%new_labels = ["ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21"]';

old_labels = ["ZT0a" "ZT3a" "ZT6a" "ZT9a" "ZT12a" "ZT15a" "ZT18a" "ZT21a" ...
              "ZT0b" "ZT3b" "ZT6b" "ZT9b" "ZT12b" "ZT15b" "ZT18b" "ZT21b"]';

new_labels = ["ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21" ... 
              "ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21"]';

tmeta = table( old_labels, new_labels);

custom_genelist = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
    "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];

custom_celltype = "tCAF";
[T1, T2] = sce_circ_phase_estimation(sce, tmeta, true, false, ...
                      [] , custom_celltype );

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

% Plot confident genes circadian expression
t0 = 0:0.1:21;
t = 0:3:21;
path = strcat(custom_celltype,'_circadian_fit_expression');
mkdir(path);
for i = 1:ngenes
    f = figure('visible','off');
    fval = T1.Amp(i).*cos( 2*pi*( t0 - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
    Rzts = table2array( T2(i,2:end) );
    plot(t0, fval);
    hold on;
    plot(t,Rzts);
    xlim([0 21])
    title("Gene - "+T1.Genes(i)+" | pvalue " + string(T1.pvalues(i)));
    xlabel('Time (hrs)');
    ylabel('Expression');
    legend({'Sine-fitted expression','Expression'},'Location','northwest');
    fname = strcat("/plot_circadian_expression_fit_", T1.Genes(i));
    fname = strcat(path,fname,".png");
    saveas(f,fname)
    hold off;
end

% Plot classic circadian genes identified with neighboars...
gidx = find(founds == 1)';
colors = colormap(jet);
path = strcat(custom_celltype,'_circadian_neighboars');
mkdir(path);
for i = gidx
    disp(i)
    icolor = 1;
    ii = 1;
    f = figure('visible','off');
    gjdx = abs( T1.Acrophase_24(:) - T1.Acrophase_24(i) ) < 0.05;
    gjdx = find( gjdx == 1)';
    disp(T2.Genes(gjdx));
    for j = gjdx
        Rzts = table2array( T2(j,2:end) );
        %disp(Rzts)
        plot(t, Rzts, 'Color', colors(icolor, :), 'Marker', 'o');
        icolor = icolor + 10;
        legend_text{ii} = T2.Genes(j);
        ii = ii + 1;
        hold on;
    end
    xlim([0 21])
    legend(legend_text, 'Location', 'bestoutside');
    title("Acrophase: "+T1.Acrophase_24(i));
    xlabel('Time (hrs)');
    ylabel('Expression');
    % Save figure
    fname = strcat("/plot_circadian_expression_neig_", T1.Genes(i));
    fname = strcat(path,fname,".png");
    saveas(f,fname)
    clear legend_text;
    % Save text file with genes
    fname = strcat("/list_circadian_expression_neig_", T1.Genes(i));
    fname = strcat(path,fname);
    tab = table(T2.Genes(gjdx));
    writetable(tab, fname,  'WriteVariableNames', 0);
    hold off;
end