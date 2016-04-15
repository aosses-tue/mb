function [SRT Stdev filenames outs] = quick_staircases(directories2check, opts, dest_folder, hoofd_folder, bSave)
% function [SRT Stdev filenames outs] = quick_staircases(directories2check, opts, dest_folder, hoofd_folder, bSave)
%
% 1. Description:
%       Since 09/03/2015 this script is compatible to extract results from
%       multiprocedure experiments.
%       Valid for 1 directory at a time
% 
% 2. Stand-alone example:
% % 2.1. Example 2015:
%       dir2check = {'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\'};
%       quick_staircases(dir2check);
% 
% % 2.2.  Example 2012-2013 (Ubuntu computer):
%       directories2check = {   'ci-Jean-Baptiste_Daumerie/20131016-LT/', ...
%                               'ci-Jean-Baptiste_Daumerie/20131022-LT/', ...
%                               'nh-Anneke_Lenssen/20130806-LT/'};
%       hoofd_folder = '~/Documenten/Meas/Meas/Experiments/Results_XML/';
%       dest_folder  = [hoofd_folder 'ci_pooled/Figures/'];
%       quick_staircases(directories2check, dest_folder, hoofd_folder);
%
% % 2.3 Example 2016 (Windows computer), using one input file:
%       file = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-4-Multiprocedure\piano-P1-P7-A4-multi-test.apr';
%       opts.bPlot = 1; % if you want to plot the figure
%       opts.bSave = 1; % if you want to save the figure
%       quick_staircases(file,opts);      
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150313_update.m
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2012-2013
% Created in    : 2012-2013
% Last update on: 30/03/2016 
% Last use on   : 30/03/2016 
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
N        = 10;

Stairs  = [];

for j = 1
    
    if isdir(directories2check)
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
        [tmp_staircase xx procID] = a3adaptiveresults(file);
        
        Nprocedures = length(fieldnames(tmp_staircase));
        for k = 1:Nprocedures
            N = max(N, length(tmp_staircase.(procID{k})));
        end
            
        % Stairs = nan(1,N); 
        Stairs = nan(Nprocedures,N);
        
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
            Stairs(1,1:length(tmp_staircase)) = tmp_staircase';
            
        else
            switch mode
                case 'mean'
                    SRT(i,1) = mean(tmp_staircase.(procID{1})(end-N4mean+1:end));
                    Stdev(i,1) = std(tmp_staircase.(procID{1})(end-N4mean+1:end));
                case 'median'
                    Md = Get_mAFC_reversals(tmp_staircase.(procID{1}));
                    SRT(i,1) = median(Md(end-N4mean+1:end));
                    Stdev(i,1) = std(Md(end-N4mean+1:end));
            end
            Stairs(1,1:length(tmp_staircase.(procID{1}))) = tmp_staircase.(procID{1})';
       
            for k = 2:Nprocedures
                switch mode
                    case 'mean'
                        exp = sprintf('SRT(i,%.0f) = mean(tmp_staircase.(procID{%.0f})(end-%.0f+1:end));',k,k,N4mean);
                        eval(exp);
                        exp = sprintf('Stdev(i,%.0f) = std(tmp_staircase.(procID{%.0f})(end-%.0f+1:end));',k,k,N4mean);
                        eval(exp);
                    case 'median'
                        exp = sprintf('Md = Get_mAFC_reversals(tmp_staircase.(procID{%}).0f);',k);
                        eval(exp);
                        exp = sprintf('SRT(i,%.0f) = median(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                        exp = sprintf('Stdev(i,%.0f) = std(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                end
                exp = sprintf('Stairs(k,1:length(tmp_staircase.(procID{%.0f}))) = tmp_staircase.(procID{%.0f});',k,k);
                eval(exp)
                
            end
        end

        if bPlot
            figure
        end
        for k = 1:Nprocedures
            if bPlot
                subplot(Nprocedures,1,k)
                plot(1:N, Stairs(k,:),'ro','LineWidth',4), hold on
                plot([0 N],[SRT(i,k) SRT(i,k)],'k--','LineWidth',2)

                plot(1:N, Stairs(k,:))
                hxLabel = xlabel('Trial number');
                set(hxLabel,'FontSize',FontSize)
                hyLabel = ylabel('Signal-to-noise ratio (dB)');
                set(hyLabel,'FontSize',FontSize)
                hLegend = legend('staircase',['SRT = ' Num2str(SRT(i,k)) 'dB; std =' Num2str(Stdev(i,k)) ' dB']);
                set(hLegend,'Location','NorthWest')
                set(hLegend,'FontSize',FontSize)

                if Nprocedures == 1
                    try
                        hTitle = title(name2figname(files(i).name));
                    catch
                        hTitle = title(name2figname(files{i}));
                    end
                else
                    try
                        tmp = fieldnames(tmp_staircase);
                        hTitle = title( [name2figname(files(i).name) ' - ' tmp{k}]);
                    catch
                        for_title = procID{k};
                        hTitle = title(for_title);
                    end
                    
                end

                try
                    set(hTitle,'FontSize',FontSize)
                end
                % ylim([ymin ymax])
                grid on
            end
            
            try
                filenames{i} = files(i).name;
            catch
                filenames{i} = files{i};
            end
            try
                name2save = strsplit(filenames{i},'.');
                name2save = name2save{1};
                name2save = strsplit(name2save,delim);
                name2save = name2save{end};
            end
        end

        if bSave
            if bPlot
                saveas(gcf,[dest_folder, name2save{1} '.eps'],'epsc');
            end
        end
        
    end
end

outs.staircases = tmp_staircase;
tmp = fieldnames(tmp_staircase);
try
    for i = 1:Nprocedures
        outs.reversals(i,:) = Get_mAFC_reversals( tmp_staircase.(tmp{i}) );
    end
end

if bSave
    disp([mfilename '.m: all files were stored at ' dest_folder])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end