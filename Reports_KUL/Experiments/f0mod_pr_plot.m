function [ACE4plot, F0m4plot] = f0mod_pr_plot(data, infodir)
% function [ACE4plot, F0m4plot] = f0mod_pr_plot(data, infodir)
%
% selection = alldata(alldata.type==LS & alldata.freq == 500,:);
% statarray = grpstats(selection,{'level'},{'mean','std'},...
%  'DataVars','response');
%
% Programmed by Tom Francart, some comments by Alejandro Osses
% Last updated: 12/06/2014 (for Windows-TU/e compatibility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(infodir,'figures_folder')
    infodir = Ensure_field(infodir,'bSave',1);
else
    infodir = Ensure_field(infodir,'bSave',0);
end

if infodir.bSave
    disp([mfilename '.m: Figures will be stored'])
else
    disp([mfilename '.m: Figures will not be stored'])
end
    
subjects    = unique(data.subject); % useful command
registers   = unique(data.register);
strategies  = unique(data.strategy);

F0m4plot = [];
ACE4plot = [];

row_ACE = find(strcmp(data.strategy,'ACE'));
dataACE = data(row_ACE,:);

row_F0m = find(strcmp(data.strategy,'F0m'));
dataF0m = data(row_F0m,:);

for is=1:length(subjects)
    
    subject = subjects{is};
    
    row_SubjectACE = find( strcmp(dataACE.subject, subject) );
    row_SubjectF0m = find( strcmp(dataF0m.subject, subject) );
    
    if length(row_SubjectACE) == length(row_SubjectF0m) && length(row_SubjectF0m) == 20;
        F0m4plot = [F0m4plot; dataF0m.percent(row_SubjectF0m)'];
        ACE4plot = [ACE4plot; dataACE.percent(row_SubjectACE)'];
    else
        disp([mfilename '.m: Subject ' subject ' appears to be incomplete...'])
        F0m4plot = [F0m4plot; nan(1,20)];
        ACE4plot = [ACE4plot; nan(1,20)];
        
        if strcmp(subject,'WD') % Completing manually Wouter's data
            ACE4plot(end,1:16)          = transpose( dataACE.percent(row_SubjectACE) );
            F0m4plot(end,[1:4, 9:16])   = transpose( dataF0m.percent(row_SubjectF0m) );
        end
    end
   
end

%% Pooled data
stPlot = label_GUI_MM; % Inline function
PlotMeans(stPlot, ACE4plot, F0m4plot); % see figPRScores_xPC

if infodir.bSave
    
    Subjects = {'S16','S14','S12','S13','S11','S15'};
    
    filename = [infodir.figures_folder 'pr-last'];
    Saveas(gcf, filename)

    for i = 1:size(ACE4plot,1)
        PlotMeans(stPlot, ACE4plot(i,:), F0m4plot(i,:));
        
        filename = [infodir.figures_folder Subjects{i} 'pr'];
        Saveas(gcf, filename)
    end
end

%% Per subject

for is=1:length(subjects)
    subject = subjects{is};
    
    figure;
    rows=3;
    cols=2;
    
    for ir=1:length(registers)
        ax=subplot(rows,cols,ir);
        hold all;
        
        for istrat=1:length(strategies)
            
            statarray = grpstats( data(strcmp(data.subject,subject) & ...
                                       data.register==registers(ir) & ...
                                       strcmp(data.strategy,strategies{istrat}) , ...
                                       : ),...
                                 {'semitones' }, ...
                                 {'mean','std'}, ...
                                  'DataVars','percent' );
            
            h=plot(statarray.semitones, statarray.mean_percent, '-x');
            
        end
        
        label_l(ax, registers(ir))%, strategies);
        
    end
    suptitle(subject);
    legend(strategies);
    
end

% %% Global (Tom's way)
% 
% figure;
% rows=3;
% cols=2;
% 
% for ir=1:length(registers)
%     subplot(rows,cols,ir);
%     hold all;
%     
%     for istrat=1:length(strategies)
%         
%         statarray = grpstats(data(...
%             data.register==registers(ir) & ...
%             strcmp(data.strategy,strategies{istrat}) ,: ),...
%             {'semitones' },{'mean','std'},...
%             'DataVars','percent');
%         
%         errorbar(statarray.semitones+(istrat-1)*0.1, statarray.mean_percent, statarray.std_percent);
%         
%     end
%     
%     label_l(ax, registers(ir))%, strategies);
%     
% end
% 
% suptitle('Grand average');
% legend(strategies);

disp([mfilename '.m: End of script'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_l(h, register) % , strategies)

plot([0 5], [50 50], ':');
set(h, 'XTick', [1:4]);
xlabel('Interval (semitones)');
ylabel('% correct');
title(register);
ylim([0 105]);
xlim([0.5 4.5]);

function paramsPlot = label_GUI_MM

try
    info = getpaths;
    addpath([info.svn.MATLAB 'Meas/Experiments/']); % get_tones_PR.m is here
    addpath([info.svn.MATLAB 'Statistics/' ]); % get_significance_level.m is here
end

% paramsPlot.SeriesLabel = {'xPC F0mod', 'xPC ACE'};
paramsPlot.SeriesLabel = {'ACE', 'F0mod'};
paramsPlot.SeriesColor = {'w'  , 'k'};
paramsPlot.scalerate = 0.9;
paramsPlot.separation_series = 0.20 * paramsPlot.scalerate;

nStrategies         = length(paramsPlot.SeriesLabel);
nRepetitions        = 9; % automate this. They were 9 presentations
nSubjects           = 6; % automate this. They were 6 subjects

paramsPlot.slev     = get_significance_level(nSubjects*nRepetitions, nStrategies); 

paramsPlot.fntsz    = ceil(10*paramsPlot.scalerate);
paramsPlot.markerSize = ceil(8*paramsPlot.scalerate);

registers           = {'G#2', 'C3', 'E3', 'G#3', 'C4'};
[tones, xTicks]     = get_tones_PR(registers);   

paramsPlot.nPlots   = length(registers);
paramsPlot.nCols    = 2;
paramsPlot.nRows    = 3;

paramsPlot.figPos   = [0 0 1024 350*paramsPlot.nRows]*paramsPlot.scalerate*0.7;

paramsPlot.xTick    = [1 2 3 4];
paramsPlot.xTickLabel = xTicks(:,2:end);
paramsPlot.yTick    = [20 40 60 80 100];
paramsPlot.xLim     = [0.5 4.5];
paramsPlot.yLim     = [0 120];

paramsPlot.XLabel   = 'Comparison notes [Hz]';
paramsPlot.YLabel   = '% correct';

for i = 1:length(registers)
    % paramsPlot.Title{i}    = [registers{i} ' (' num2str(xTicks(i,1)) ' - ' num2str(xTicks(i,end)) 'Hz)']; % Matthias' format
    paramsPlot.Title{i}    = [registers{i} ' (' num2str(xTicks(i,1)) ' Hz)'];
end

