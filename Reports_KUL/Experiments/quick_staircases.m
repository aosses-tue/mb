function quick_staircases(directories2check, dest_folder, hoofd_folder)
% function quick_staircases(directories2check, dest_folder, hoofd_folder)
%
% % Example:
%   directories2check = {   'ci-Jean-Baptiste_Daumerie/20131016-LT/', ...
%                           'ci-Jean-Baptiste_Daumerie/20131022-LT/', ...
%                           'nh-Anneke_Lenssen/20130806-LT/'};
%   hoofd_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
%   dest_folder  = [hoofd_folder 'ci_pooled/Figures/'];
%   quick_staircases(directories2check, dest_folder, hoofd_folder);
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bSave = 1;

FontSize = 14;
N       = 20;
Stairs  = [];
ymin    = -5;
ymax    = 15;

for j = 1:length(directories2check)
    
    files = dir([hoofd_folder directories2check{j} '*.apr']);
    for i = 1:length(files)
        try
            file = [hoofd_folder directories2check{j} files(i).name];
            tmp_staircase = a3adaptiveresults(file);
            SRT = mean(tmp_staircase(end-6+1:end));

            Stairs = [Stairs; nan(1,N)];
            Stairs(end,1:length(tmp_staircase)) = tmp_staircase';

            figure
            plot(1:N, Stairs(end,:),'ro','LineWidth',4), hold on
            plot([0 N],[SRT SRT],'k--','LineWidth',2)
            plot(1:N, Stairs(end,:))
            hxLabel = xlabel('Trial number');
            set(hxLabel,'FontSize',FontSize)
            hyLabel = ylabel('Signal-to-noise ratio (dB)');
            set(hyLabel,'FontSize',FontSize)
            hLegend = legend('staircase',['SRT = ' Num2Str(SRT) 'dB; std =' Num2Str(std(tmp_staircase(end-6+1:end))) ' dB']);
            set(hLegend,'Location','NorthWest')
            set(hLegend,'FontSize',FontSize)

            hTitle = title(name2figname(files(i).name));
            set(hTitle,'FontSize',FontSize)
            ylim([ymin ymax])
            grid on

            name2save = strsplit(files(i).name,'.');

            if bSave
                %saveas(gcf,[dest_folder, name2save{1} '.eps'],'epsc');
            end
        catch
            warning(['File ' file ' seems not to correspond to the result of an adaptive procedure... '])
        end
        
    end
end

if bSave
    disp([mfilename '.m: all files were stored at ' dest_folder])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end