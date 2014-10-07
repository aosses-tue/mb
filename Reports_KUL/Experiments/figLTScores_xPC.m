function [h p] = figLTScores_xPC(file, bPlot, stPlot)
% function [h p] = figLTScores_xPC(file, bPlot, stPlot)
%
% Processing data for Sentences-in-noise (LIST) conducted on:
%       - NH Subjects (July-August 2013) = NH-Sentence-recognition.txt
%       - NH Subjects (30 August 2013)   = NH-Sentence-recognition-retest.txt
%       - CI Subject: Romain (19-Sept-2013)
%
% file  -   by default: NH-Sentence-recognition.txt (data presented in 
%           Predoctoral thesis)
% p.aceData - Results for ACE strategy (strategy = 0)
% p.f0mData - Results for F0mod strategy (strategy = 1)
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List_of_files = {   'confidenceInterval.m', ...
%                     'get_significance_level.m', ...
%                     'setupPlotConf.m'};
%                 
% [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
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
STD_dB = 2.5;

import_scores_from_txt(filename);
data        = evalin('base','data');
textdata    = evalin('base','textdata');

if strcmp(file, 'NH-Sentence-recognition.txt')
    rows2delete = 23:26; % Jonas' data
    data(rows2delete,:) = [];
end

column_strategy = 1; 
column_num_list = 3;
column_STD      = 4;
colun_subject   = 5;

try
    rows2delete     = find(data(:,column_STD)>STD_dB);
    disp([mfilename '.m: excluding data with STD greater than ' num2str(STD_dB) ' dB'])
    pause(1)
    data(rows2delete,2:4) = nan;
end

idx_ACEown      = find(data(:,column_strategy)==10); 
idx_ACE         = find(data(:,column_strategy)==0);
idx_F0m         = find(data(:,column_strategy)==1);
idx_F0m_FIR_EE  = find(data(:,column_strategy)==2);

column_score_sentence   = 2; % after deleting

aceDataOwn      = data(idx_ACEown     , [column_score_sentence]); 
aceData         = data(idx_ACE        , [column_score_sentence]); 
f0mData         = data(idx_F0m        , [column_score_sentence]);
f0mFIREEData    = data(idx_F0m_FIR_EE , [column_score_sentence]);

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
    stPlot              = Ensure_field(stPlot,'YLim',[-10 10]);
    stPlot.yTick        = min(stPlot.YLim)+2:2:max(stPlot.YLim)-2;
    stPlot              = Ensure_field(stPlot,'xTickLabel',{''});
    stPlot.YGrid        = 'on';
    
    stPlot.TitleSuffix  = {''};
    stPlot              = Ensure_field(stPlot,'SeriesLabel',{'xPC F0mod', 'xPC ACE', 'ACE baseline'});
    stPlot              = Ensure_field(stPlot,'TitleHead','LIST material'); 
    stPlot.YLabel       = 'SNR (dB)';
    stPlot              = Ensure_field(stPlot,'XLabel',' ');
    stPlot.ReverseData  = 1; % from worst to best score
    stPlot.separation_series = 0.1;
    
    if length(aceDataOwn)~=0
        PlotMeans3series(stPlot, f0mData, aceData, aceDataOwn);
    else
        SeriesLabel         = {'xPC F0mod', 'xPC ACE'};
        PlotMeans(stPlot, f0mData, aceData);
    end

    h = gcf;
    set(gcf, 'PaperPositionMode','auto')
end

p.aceData = data(idx_ACE,:);
p.f0mData = data(idx_F0m,:);
p.column_subject    = 5;
p.column_score      = 2;
p.column_strategy   = 1;
try
    p.f0FIREEData = data(idx_F0m_FIR_EE,:);
end

try
    p.aceDataOwn = data(idx_ACEown,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end