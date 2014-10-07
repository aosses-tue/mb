function experiment_report_20131030_JBD_staircases
% function experiment_report_20131030_JBD_staircases
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
dest_folder     = [hoofd_folder 'ci_pooled/Figures/'];

SubjectName = 'CIs';

bPlot = 1;
bSave = 1;

directory = {'ci-Jean-Baptiste_Daumerie/20131016-LT/', ...
             'ci-Jean-Baptiste_Daumerie/20131022-LT/'};

N       = 25;
Stairs  = [];
ymin    = -5;
ymax    = 15;

for j = 1:length(directory)
    files = dir([hoofd_folder directory{j} '*.apr']);
    for i = 1:length(files)
        tmp_staircase = a3adaptiveresults([hoofd_folder directory{j} files(i).name]);
        SRT = mean(tmp_staircase(end-6+1:end));
        
        Stairs = [Stairs; nan(1,N)];
        Stairs(end,1:length(tmp_staircase)) = tmp_staircase';
        
        figure
        plot(1:25, Stairs(end,:),'ro','LineWidth',4), hold on
        plot([0 25],[SRT SRT],'k--','LineWidth',2)
        plot(1:25, Stairs(end,:))
        xlabel('Trial number')
        ylabel('Signal-to-noise ratio (dB)')
        legend('staircase',['SRT = ' Num2Str(SRT) 'dB; std =' Num2Str(std(tmp_staircase(end-6+1:end))) ' dB'])
        title(name2figname(files(i).name))
        ylim([ymin ymax])
        grid on
        
        saveas(gcf,[dest_folder,files(i).name '.eps'],'epsc');
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end