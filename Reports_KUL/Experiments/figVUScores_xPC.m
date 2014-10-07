function [h p] = figVUScores_xPC(file, bPlot, stPlot)
% function [h p] = figVUScores_xPC(file, bPlot, stPlot)
%
% Processing data for Sentences-in-noise (VU) conducted on:
%       - CI Subjects (19 September 2013)   = YYYYYYYYYYYYYYY.txt
%
% file  -   by default: YYYYYYYYYYYYYYYYY.txt (data presented in 
%           Predoctoral thesis)
% p.aceData - Results for ACE strategy (strategy = 0)
% p.f0mData - Results for F0mod strategy (strategy = 1)
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'confidenceInterval.m', ...
                    'get_significance_level.m', ...
                    'setupPlotConf.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    file = 'NH-Sentence-recognition.txt';
    filepath    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
    filename    = [filepath file];
else
    filename = file;
end

if nargin < 2
    bPlot = input([mfilename '.m - Do you want to generate the plots for the LIST scores? (1 = yes / 0 = no): ']);
end

h = []; % initialisation of handle

import_scores_from_txt(filename);
data        = evalin('base','data');
textdata    = evalin('base','textdata');

column_strategy = 1; 
column_num_list = 3;

idx_ACEown      = find(data(:,column_strategy)==10); 
idx_ACE         = find(data(:,column_strategy)==0);
idx_F0m         = find(data(:,column_strategy)==1);
idx_F0m_FIR_EE  = find(data(:,column_strategy)==2);

columns2discard = [column_strategy column_num_list];

column_score_sentence   = 2; % after deleting

aceDataOwn      = data(idx_ACEown     , column_score_sentence); 
aceData         = data(idx_ACE        , column_score_sentence); 
f0mData         = data(idx_F0m        , column_score_sentence);
f0mFIREEData    = data(idx_F0m_FIR_EE , column_score_sentence);

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPlot == 1
    stPlot.nCols        = 1;
    stPlot.nRows        = 1;
    stPlot.xGap         = 15;
    stPlot.yGap         = 80;
    stPlot.fntsz        = 16;
    stPlot.markerSize   = 12;
    stPlot.margins      = [70 70 50 90];
    stPlot.figPos       = [0 0 1024 450]; % 225 each row % left, right, top, bottom
    stPlot.xLim         = [  0  2];
    stPlot = Ensure_field(stPlot,'yLim', [0 20]);
    stPlot.yTick        = stPlot.yLim(1)+2:2:stPlot.yLim(2)-2;
    stPlot.xTickLabel   = {''};
    stPlot.YGrid        = 'on';
    % 
    
    stPlot.TitleSuffix  = {''};
    stPlot.TitleHead    = ['VU material']; 
    stPlot.YLabel       = 'SNR (dB)';
    stPlot.XLabel       = ' ';
    stPlot.ReverseData  = 1; % from worst to best score
    stPlot.separation_series = 0.1;

    try
        stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod', 'ACE baseline'};
        PlotMeans3series(stPlot, aceData, f0mData,  aceDataOwn);
    catch
        stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};
        PlotMeans(stPlot, aceData, f0mData);
    end
    
    h = gcf;
    set(gcf, 'PaperPositionMode','auto')
end

p.aceData = data(idx_ACE,:);
p.f0mData = data(idx_F0m,:);
try
    p.f0FIREEData = data(idx_F0m_FIR_EE,:);
end

try
    p.aceDataOwn = data(idx_ACEown,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end