
    %Grid for more resolution
    t = 0:0.5:time_cycle;
    fval = fmdl.amp * cos( 2*pi*( t - fmdl.acro)/24 ) + fmdl.mesor; % Period 24
    %plotit = true;
    if plotit
        figure;
        plot(t,fval);
        xticks(0:3:21);    
        xticklabels(strrep(cL(1:end),'a',''));
        xlim([-0.1 21.1])
        title(gs,cellt);
        hold on;
        plot(time_grid,R,'-s');
        hold on;
        % Shadow
        a=ylim;
        ylim([a(1) a(2)*1.2]);
        a=ylim;
        x=[10.5 21 21 10.5];
        y=[0 0 a(2) a(2)];
        patch('Faces',[1 2 3 4],'Vertices',[x' y'], ...
        'FaceColor','black','FaceAlpha',.3,'Edgecolor','none');
        legend('sine-fit','Data',"Night")
        hold on;
        saveas(gcf,strcat(gs,cellt),'png')
        close
    end