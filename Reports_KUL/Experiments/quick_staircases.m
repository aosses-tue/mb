function SRT = quick_staircases(directories2check, dest_folder, hoofd_folder, bSave)
% function SRT = quick_staircases(directories2check, dest_folder, hoofd_folder, bSave)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments
% 
% 2. Stand-alone example:
% % Example 2015:
%   dir2check = {'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\'};
%   quick_staircases(dir2check);
% 
% % Example 2012-2013 (Ubuntu computer):
%   directories2check = {   'ci-Jean-Baptiste_Daumerie/20131016-LT/', ...
%                           'ci-Jean-Baptiste_Daumerie/20131022-LT/', ...
%                           'nh-Anneke_Lenssen/20130806-LT/'};
%   hoofd_folder = '~/Documenten/Meas/Meas/Experiments/Results_XML/';
%   dest_folder  = [hoofd_folder 'ci_pooled/Figures/'];
%   quick_staircases(directories2check, dest_folder, hoofd_folder);
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2012-2013
% Created in    : 2012-2013
% Last update on: 09/03/2015 % Update this date manually
% Last use on   : 09/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    bSave = 0;
end

if nargin < 3
    hoofd_folder = '';
end

if nargin < 2
    dest_folder = Get_TUe_paths('outputs');
end

SRT = [];

FontSize = 14;
N       = 20;
Stairs  = [];
ymin    = -5;
ymax    = 15;

for j = 1:length(directories2check)
    
    files = dir([hoofd_folder directories2check{j} '*.apr*']);
    for i = 1:length(files)
        try
            file = [hoofd_folder directories2check{j} files(i).name];
            tmp_staircase = a3adaptiveresults(file);
            Stairs = [Stairs; nan(1,N)];
            NumProcedures = 1;
            
            if ~isstruct( tmp_staircase );
                SRT = mean(tmp_staircase(end-6+1:end));
                Stairs(end,1:length(tmp_staircase)) = tmp_staircase';
                Stdev = std(tmp_staircase(end-6+1:end));
       
            else
                SRT = mean(tmp_staircase.procedure1(end-6+1:end));
                Stairs(end,1:length(tmp_staircase.procedure1)) = tmp_staircase.procedure1';
                Stdev = std(tmp_staircase.procedure1(end-6+1:end));
                
                for k = 2:length(fieldnames(tmp_staircase)) 
                    exp = sprintf('SRT(%.0f) = mean(tmp_staircase.procedure%.0f(end-6+1:end));',k,k);
                    eval(exp);
                    exp = sprintf('Stairs(end,1:length(tmp_staircase.procedure%.0f)) = tmp_staircase.procedure%.0f;',k,k);
                    eval(exp)
                    exp = sprintf('Stdev(%.0f) = std(tmp_staircase.procedure%.0f(end-6+1:end));',k,k);
                    eval(exp);
                    NumProcedures = NumProcedures + 1;
                end
            end

            figure
            for k = 1:NumProcedures
                subplot(NumProcedures,1,k)
                plot(1:N, Stairs(end,:),'ro','LineWidth',4), hold on
                
                plot([0 N],[SRT(k) SRT(k)],'k--','LineWidth',2)
                
                plot(1:N, Stairs(end,:))
                hxLabel = xlabel('Trial number');
                set(hxLabel,'FontSize',FontSize)
                hyLabel = ylabel('Signal-to-noise ratio (dB)');
                set(hyLabel,'FontSize',FontSize)
                hLegend = legend('staircase',['SRT = ' Num2str(SRT(k)) 'dB; std =' Num2str(Stdev(k)) ' dB']);
                set(hLegend,'Location','NorthWest')
                set(hLegend,'FontSize',FontSize)

                if NumProcedures == 1
                    hTitle = title(name2figname(files(i).name));
                else
                    tmp = fieldnames(tmp_staircase);
                    hTitle = title( [name2figname(files(i).name) ' - ' tmp{k}]);
                end
                
                set(hTitle,'FontSize',FontSize)
                ylim([ymin ymax])
                grid on

                name2save = strsplit(files(i).name,'.');
            end

            if bSave
                saveas(gcf,[dest_folder, name2save{1} '.eps'],'epsc');
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