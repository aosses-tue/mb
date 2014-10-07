function h = figMCScores_xPC(dirExpRes, instr)
% h = function figMCScores_xPC(dirExpRes, instr)
% 
% Get plots for MCI Experiments. 
% MCI plot (global result) is similar to Figure  3.10, page 64 (but different format)
% Matthias' thesis. 
%
% Examples:
%       figMCScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'UW'})
%
% Programmed by Alejandro, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'confidenceInterval.m', ...
                    'get_confusion_matrix.m', ...
                    'get_mci_contours.m', ...
                    'plot_MCI_results.m', ...
                    'pool_mci_results.m', ...
                    'process_stimuli_exp.m', ...
                    'setupPlotConf.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dirExpRes   = '/home/alejandro/Documenten/MM/data/Experiments/Results_XML_Strategies/';
end

if nargin < 2
    instr{1} = '';
end

% Each row of aceMCIData corresponds to a different subject
% Columns from 1:9, 10:18, 19:27, 28:36 represent the scores for each sequence
% considering a separation between tones of 1, 2, 3, 4 semitones respectively
% [aceMCIData, f0mMCIData] = readMCIData; % Matthias' format

regs        = {'G#2'};

[aceMCIData, f0mMCIData, p] = readMCData(dirExpRes, regs, instr{1});
[nSubjects nIntervalsTotal] = size(aceMCIData);
display(['Number of subjects taken into account: ', num2str(nSubjects)])
subjects = p.subjects;

nRoving     = 5; % The 9 contours are being presented 5 times
nChoices    = 9; % Number of contours
nIntervals  = nIntervalsTotal/nChoices;

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Number of subjects plotted: ' num2str(nSubjects)])
disp('Subjects/Experiment folder:  ')
for count = 1:nSubjects
    disp([subjects{count,1} '/' subjects{count,2}]);
end

meanACE     = mean(aceMCIData);
yDataACE    = reshape(meanACE, nChoices, nIntervals); % 9x4, each row is a different contour averaged across subjects

meanF0m     = mean(f0mMCIData);
yDataF0m    = reshape(meanF0m, nChoices, nIntervals); 

nTotalTrials= nSubjects*nChoices*nRoving;
stPlot.slev = get_significance_level(nTotalTrials, nChoices); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margins     = [50 50 25 50];
figPos      = [50 50 1024 768/2];
stPlot.xLim = [0.8 4.2]; % Manual setup
stPlot.yLim = [-30 120]; % -10
stPlot.xTick= 1:4; % Defines the x-axis on plot, size: 1x4 
stPlot.yTick= 20:20:100;
doubleAxis  = 0;
nRows       = 1;
nCols       = 1;
xGap        = 0;
yGap        = 0;
offset_ACE  = -0.08;
offset_f0m  = -offset_ACE; % offset around xTick

stPlot.SeriesLabel  = {'ACE','F0mod'};
stPlot.XLabel       = 'Interval (semitones)';
stPlot.YLabel       = '% Correct';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stPlot.xTick    = 1:4;
stPlot.LocationLegend = 'NorthWest';
stPlot.TitleHead = 'MCI Results';
PlotMeans(yDataACE, yDataF0m, stPlot);
h(1) = gcf;

stPlot.figPos   = [50 50 1024 768];
stPlot.xTick    = 1:9;
stPlot.xLim     = [0.8 10]; % Manual setup
stPlot.nCols    = 2;
stPlot.nRows    = 2;
stPlot.TitleHead = '';
stPlot.nPlots   = stPlot.nCols*stPlot.nRows;
stPlot.XLabel       = '';
stPlot.xTickLabel = get_mci_contours;

for i=1:stPlot.nPlots
    stPlot.TitleSuffix{1,i} = ['Interval group ' int2str(i)];
end
PlotMeans(aceMCIData, f0mMCIData, stPlot);
h(2) = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end

end
