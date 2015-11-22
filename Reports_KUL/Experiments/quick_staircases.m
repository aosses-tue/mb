function [SRT Stdev filenames] = quick_staircases(directories2check, opts, dest_folder, hoofd_folder, bSave)
% function [SRT Stdev filenames] = quick_staircases(directories2check, opts, dest_folder, hoofd_folder, bSave)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments.
%       Valid for 1 directory at a time
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
% Last update on: 13/03/2015 
% Last use on   : 13/03/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    bSave = 0;
end

if nargin < 4
    hoofd_folder = '';
end

if nargin < 3
    dest_folder = Get_TUe_paths('outputs');
end

if nargin < 2
    opts = [];
end

opts = Ensure_field(opts,'mode','mean'); % 'median'
if nargout == 0
    opts = Ensure_field(opts,'bPlot',1);
else
    opts = Ensure_field(opts,'bPlot',0);
end
opts = Ensure_field(opts,'N4mean',6);
opts = Ensure_field(opts,'filter','');

mode    = opts.mode;
N4mean  = opts.N4mean;
bPlot   = opts.bPlot;

SRT = [];

FontSize = 14;
N       = 25;

Stairs  = [];
% ymin    = -5;
% ymax    = 15;

for j = 1
    
    if isdir(directories2check)
        % files = dir([hoofd_folder directories2check{j} '*.apr*']);
        files = dir([hoofd_folder directories2check opts.filter '*.apr*']);
    else
        files{1} = directories2check;
    end
    
    for i = 1:length(files)
        
        if isdir(directories2check)
            % file = [hoofd_folder directories2check{j} files(i).name];
            file = [hoofd_folder directories2check files(i).name];
        else
            file = files{1};
        end
        tmp_staircase = a3adaptiveresults(file);
        Stairs = [Stairs; nan(1,N)];
        NumProcedures = 1;

        if ~isstruct( tmp_staircase );
            
            switch mode
                case 'mean'
                    SRT(i,1) = mean(tmp_staircase(end-N4mean+1:end));
                    Stdev = std(tmp_staircase(end-N4mean+1:end));
                case 'median'
                    Md = Get_mAFC_reversals(tmp_staircase);
                    SRT(i,1) = median(Md);
                    Stdev = std(tmp_staircase(end-N4mean+1:end));
            end
            Stairs(end,1:length(tmp_staircase)) = tmp_staircase';
            
        else
            switch mode
                case 'mean'
                    SRT(i,1) = mean(tmp_staircase.procedure1(end-N4mean+1:end));
                    Stdev(i,1) = std(tmp_staircase.procedure1(end-N4mean+1:end));
                case 'median'
                    Md = Get_mAFC_reversals(tmp_staircase.procedure1);
                    SRT(i,1) = median(Md(end-N4mean+1:end));
                    Stdev(i,1) = std(Md(end-N4mean+1:end));
            end
            Stairs(end,1:length(tmp_staircase.procedure1)) = tmp_staircase.procedure1';
            

            for k = 2:length(fieldnames(tmp_staircase)) 
                switch mode
                    case 'mean'
                        exp = sprintf('SRT(i,%.0f) = mean(tmp_staircase.procedure%.0f(end-%.0f+1:end));',k,k,N4mean);
                        eval(exp);
                        exp = sprintf('Stdev(i,%.0f) = std(tmp_staircase.procedure%.0f(end-%.0f+1:end));',k,k,N4mean);
                        eval(exp);
                    case 'median'
                        exp = sprintf('Md = Get_mAFC_reversals(tmp_staircase.procedure%.0f);',k);
                        eval(exp);
                        exp = sprintf('SRT(i,%.0f) = median(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                        exp = sprintf('Stdev(i,%.0f) = std(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                end
                exp = sprintf('Stairs(end,1:length(tmp_staircase.procedure%.0f)) = tmp_staircase.procedure%.0f;',k,k);
                eval(exp)
                
                NumProcedures = NumProcedures + 1;
            end
        end

        if bPlot
            figure
        end
        for k = 1:NumProcedures
            if bPlot
                subplot(NumProcedures,1,k)
                plot(1:N, Stairs(end,:),'ro','LineWidth',4), hold on
                plot([0 N],[SRT(i,k) SRT(i,k)],'k--','LineWidth',2)

                plot(1:N, Stairs(end,:))
                hxLabel = xlabel('Trial number');
                set(hxLabel,'FontSize',FontSize)
                hyLabel = ylabel('Signal-to-noise ratio (dB)');
                set(hyLabel,'FontSize',FontSize)
                hLegend = legend('staircase',['SRT = ' Num2str(SRT(i,k)) 'dB; std =' Num2str(Stdev(i,k)) ' dB']);
                set(hLegend,'Location','NorthWest')
                set(hLegend,'FontSize',FontSize)

                if NumProcedures == 1
                    hTitle = title(name2figname(files(i).name));
                else
                    tmp = fieldnames(tmp_staircase);
                    hTitle = title( [name2figname(files(i).name) ' - ' tmp{k}]);
                end

                set(hTitle,'FontSize',FontSize)
                % ylim([ymin ymax])
                grid on
            end
            filenames{i} = files(i).name;
            name2save = strsplit(files(i).name,'.');
        end

        if bSave
            if bPlot
                saveas(gcf,[dest_folder, name2save{1} '.eps'],'epsc');
            end
        end
        
    end
end

if bSave
    disp([mfilename '.m: all files were stored at ' dest_folder])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end