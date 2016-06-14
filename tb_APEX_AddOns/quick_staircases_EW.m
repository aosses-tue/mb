function [AdaptiveParam Stdev filenames] = quick_staircases(file, opts, dest_folder)
% function [AdaptiveParam Stdev filenames] = quick_staircases(file, opts, dest_folder)
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
% Last update on: 18/04/2016 
% Last use on   : 18/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    bSave = 0;
else
    bSave = 1;
end

if nargin < 2
    opts = [];
end

opts = ef(opts,'mode','mean'); % 'median'
if nargout == 0
    opts = ef(opts,'bPlot',1);
else
    opts = ef(opts,'bPlot',0);
end
opts = ef(opts,'N4mean',4);
opts = ef(opts,'filter','');

mode    = opts.mode;
N4mean  = opts.N4mean;
bPlot   = opts.bPlot;

AdaptiveParam = [];

FontSize = 14;
N        = 10;

Stairs  = [];

for j = 1
    
    for i = 1 % only one file
        
        % file = files{1};
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
                    try
                        AdaptiveParam(i,1) = mean(tmp_staircase(end-N4mean+1:end));
                        Stdev = std(tmp_staircase(end-N4mean+1:end));
                    catch
                        AdaptiveParam(i,1) = nan;
                        Stdev = nan;
                    end
                    
                case 'median'
                    Md = Get_mAFC_reversals(tmp_staircase);
                    try
                        AdaptiveParam(i,1) = median(Md);
                        Stdev = std(tmp_staircase(end-N4mean+1:end));
                    catch
                        AdaptiveParam(i,1) = nan;
                        Stdev = nan;
                    end
                    
            end
            Stairs(1,1:length(tmp_staircase)) = tmp_staircase';
            
        else
            switch mode
                case 'mean'
                    try
                        AdaptiveParam(i,1) = mean(tmp_staircase.(procID{1})(end-N4mean+1:end));
                        Stdev(i,1) = std(tmp_staircase.(procID{1})(end-N4mean+1:end));
                    catch
                        AdaptiveParam(i,1) = nan;
                        Stdev(i,1) = nan;
                    end
                    
                case 'median'
                    Md = Get_mAFC_reversals(tmp_staircase.(procID{1}));
                    AdaptiveParam(i,1) = median(Md(end-N4mean+1:end));
                    Stdev(i,1) = std(Md(end-N4mean+1:end));
            end
            Stairs(1,1:length(tmp_staircase.(procID{1}))) = tmp_staircase.(procID{1})';
       
            for k = 2:Nprocedures
                switch mode
                    case 'mean'
                        try
                            exp = sprintf('AdaptiveParam(i,%.0f) = mean(tmp_staircase.(procID{%.0f})(end-%.0f+1:end));',k,k,N4mean);
                            eval(exp);
                            exp = sprintf('Stdev(i,%.0f) = std(tmp_staircase.(procID{%.0f})(end-%.0f+1:end));',k,k,N4mean);
                            eval(exp);
                        catch
                            exp = sprintf('AdaptiveParam(i,%.0f) = nan;',k);
                            eval(exp);
                            exp = sprintf('Stdev(i,%.0f) = nan;',k);
                            eval(exp);
                        end
                            
                    case 'median'
                        exp = sprintf('Md = Get_mAFC_reversals(tmp_staircase.(procID{%}).0f);',k);
                        eval(exp);
                        exp = sprintf('AdaptiveParam(i,%.0f) = median(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                        exp = sprintf('Stdev(i,%.0f) = std(Md(end-%.0f+1:end));',k,N4mean);
                        eval(exp);
                end
                exp = sprintf('Stairs(k,1:length(tmp_staircase.(procID{%.0f}))) = tmp_staircase.(procID{%.0f});',k,k);
                eval(exp)
                
            end
        end

        for k = 1:Nprocedures
            if bPlot
                
                if ~isnan(AdaptiveParam(i,k))
                    figure;
                    plot(1:N, Stairs(k,:),'ro','LineWidth',4), hold on
                    plot([0 N],[AdaptiveParam(i,k) AdaptiveParam(i,k)],'k--','LineWidth',2)

                    plot(1:N, Stairs(k,:))
                    hxLabel = xlabel('Trial number');
                    set(hxLabel,'FontSize',FontSize)
                    hyLabel = ylabel('Adaptive parameter');
                    set(hyLabel,'FontSize',FontSize)
                    hLegend = legend('staircase',['AdaptiveParam = ' num2str(AdaptiveParam(i,k)) 'dB; std =' num2str(Stdev(i,k)) ' dB']);
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
            end
            
            try
                name2save = strsplit(file,'.');
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
        
        if nargin > 0
            
            for k = 1:Nprocedures
                
                fprintf( '\t Adaptive parameter = %.4f \t (proc.: %s, \t avg of last %.0f trials) \n',AdaptiveParam(k),procID{k},opts.N4mean);
                
            end
            
        end
        
    end
end

if bSave
    disp([mfilename '.m: all files were stored at ' dest_folder])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end