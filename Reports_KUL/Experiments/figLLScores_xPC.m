function [h p] = figLLScores_xPC(filename, stPlot)
% function [h p] = figLLScores_xPC(filename, stPlot)
%
% Strategy = 0 = ACE
% Strategy = 1 = F0mod
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List_of_files = {   'confidenceInterval.m', ...
%                     'get_significance_level.m', ...
%                     'setupPlotConf.m'};
%                 
% [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    file        = 'NH-Word-recognition.txt';
    directory   = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
    filename = [directory file];
    column_subject          = 1;
    column_strategy         = 2; 
    column_score_word       = 3;
    column_score_phoneme    = 5;
    column_score_cons       = 7;
    column_score_vowel      = 6;
    column_num_list         = 4;
else
    column_subject          = 1;
    column_SNR              = 2;
    column_strategy         = 3; 
    column_score_word       = column_strategy + 1;
    column_score_phoneme    = column_strategy + 3;
    column_score_cons       = column_strategy + 5;
    column_score_vowel      = column_strategy + 4;
    column_num_list         = column_strategy + 2;  
    column_report           = column_strategy + 6;
end

if ~exist('stPlot','var')
    stPlot = [];
end

import_scores_from_txt(filename);
data = evalin('base','data');
textdata = evalin('base','textdata');

idx_ACEown = find(data(:,column_strategy)==10); % baseline
idx_ACE = find(data(:,column_strategy)==0);
idx_F0m = find(data(:,column_strategy)==1);

aceDataOwn = data(idx_ACEown,[column_score_word column_score_phoneme])*100;
aceData = data(idx_ACE,[column_score_word column_score_phoneme])*100; % in percentage. out of: 10 words
f0mData = data(idx_F0m,[column_score_word column_score_phoneme])*100;

column_scores   = 1; % after deleting
column_scores_per_phoneme = 2; % after deleting

aceDataOwn(:,column_scores)           = aceDataOwn(:,column_scores)/10;
aceDataOwn(:,column_scores_per_phoneme) = aceDataOwn(:,column_scores_per_phoneme)/30;
aceData(:,column_scores)              = aceData(:,column_scores)/10;
aceData(:,column_scores_per_phoneme)  = aceData(:,column_scores_per_phoneme)/30;
f0mData(:,column_scores)              = f0mData(:,column_scores)/10;
f0mData(:,column_scores_per_phoneme)  = f0mData(:,column_scores_per_phoneme)/30;

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.nStrategies       = 2;
p.nSemitones        = 1;
p.nReferences       = 2;
nStrategies         = p.nStrategies;
nSemitones          = p.nSemitones;

if length(stPlot) ~= 0

    stPlot.nPlots       = p.nReferences;
    stPlot.nCols        = 2;
    stPlot.nRows        = 1;
    stPlot.xGap         = 15;
    stPlot.yGap         = 80;
    stPlot.fntsz        = 16;
    stPlot.markerSize   = 12;
    stPlot.margins      = [70 70 50 90];
    stPlot.figPos       = [0 0 1024 450];   % 225 each row
                                            % left, right, top, bottom

    stPlot      = Ensure_field(stPlot,'xTick',1:nSemitones);
    stPlot      = Ensure_field(stPlot,'yTick',0:20:100);
    stPlot.xTickLabel = {''};
    stPlot.YGrid = 'on';

    SeriesLabel         = {'ACE', 'F0mod'};
    stPlot.TitleSuffix  = {' scores per word', ' scores per phoneme'};
    stPlot.SeriesLabel  = SeriesLabel;
    stPlot.Title        = {'Lilliput list: '}; 
    stPlot.YLabel       = '% Correct';
    stPlot.XLabel       = '';
    stPlot.ReverseData = 1; % from worst to best score
    stPlot.xLim = [0 2];
    stPlot.separation_series = 0.1;

    PlotMeans(stPlot, aceData, f0mData);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = gcf;
    set(gcf, 'PaperPositionMode','auto')
    
else
    h = [];
end

if ~exist('column_SNR','var')
    p.aceData = data(idx_ACE, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_subject column_report]);
    p.f0mData = data(idx_F0m, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_subject column_report]);
    if length(idx_ACEown) ~= 0
        p.aceDataOwn = data(idx_ACEown, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_subject column_report]);
    end
    p.column_subject = 5;
else
    p.aceData = data(idx_ACE, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_SNR column_subject column_report]);
    p.f0mData = data(idx_F0m, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_SNR column_subject column_report]);
    if length(idx_ACEown) ~= 0
        p.aceDataOwn = data(idx_ACEown, [column_score_word column_score_phoneme column_score_cons column_score_vowel column_SNR column_subject column_report]);
    end
    p.column_SNR = 5;
    p.column_subject = 6;
end
p.column_report         = p.column_subject + 1;
p.column_score_word     = 1;
p.column_score_phoneme  = 2;
p.column_score_cons     = 3;
p.column_score_vowel    = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
