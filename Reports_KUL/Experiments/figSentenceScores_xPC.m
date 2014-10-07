function [h p] = figSentenceScores_xPC(filename, bPlot)
% function [h p] = figSentenceScores_xPC(filename, bPlot)
%
% Sentences at fixed SNR
%
% Strategy = 0 = ACE
% Strategy = 1 = F0mod
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'confidenceInterval.m', ...
                    'get_significance_level.m', ...
                    'setupPlotConf.m'};

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

column_SNR              = 2;
column_strategy         = 3; 
column_score_sentence   = column_strategy + 1;
column_score_word       = column_strategy + 3;
column_num_list         = column_strategy + 2;   

if nargin < 2
    bPlot = 0;
end
import_scores_from_txt(filename);
data = evalin('base','data');
textdata = evalin('base','textdata');

idx_ACEown = find(data(:,column_strategy)==10); % baseline
idx_ACE = find(data(:,column_strategy)==0);
idx_F0m = find(data(:,column_strategy)==1);

aceDataOwn = data(idx_ACEown,[column_score_sentence column_score_word])*100;
aceData = data(idx_ACE,[column_score_sentence column_score_word])*100; % in percentage. out of: 10 words
f0mData = data(idx_F0m,[column_score_sentence column_score_word])*100;

column_scores_sentence_new  = 1; % after deleting
column_scores_word_new      = 2; % after deleting


aceDataOwn(:,column_scores_sentence_new)           = aceDataOwn(:,column_scores_sentence_new)/10;
aceDataOwn(:,column_scores_word_new) = aceDataOwn(:,column_scores_word_new);
aceData(:,column_scores_sentence_new)              = aceData(:,column_scores_sentence_new)/10;
aceData(:,column_scores_word_new)  = aceData(:,column_scores_word_new);
f0mData(:,column_scores_sentence_new)              = f0mData(:,column_scores_sentence_new)/10;
f0mData(:,column_scores_word_new)  = f0mData(:,column_scores_word_new);

if bPlot == 1
    % Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p.nStrategies       = 2;
    p.nSemitones        = 1;
    p.nReferences       = 2;
    nStrategies         = p.nStrategies;
    nSemitones          = p.nSemitones;
    nReferences         = p.nReferences;        stPlot.nPlots       = nReferences;

    stPlot.nCols        = 2;
    stPlot.nRows        = 1;
    stPlot.xGap         = 15;
    stPlot.yGap         = 80;
    stPlot.fntsz        = 16;
    stPlot.markerSize   = 12;
    stPlot.margins      = [70 70 50 90];
    stPlot.figPos       = [0 0 1024 450];   % 225 each row
                                            % left, right, top, bottom

    stPlot.xTick = 1:nSemitones;
    stPlot.yTick = 0:20:100;
    stPlot.xTickLabel = {''};
    stPlot.YGrid = 'on';

    SeriesLabel = {'ACE', 'F0mod'};
    stPlot.TitleSuffix = {' scores per word', ' scores per phoneme'};
    stPlot.SeriesLabel = SeriesLabel;
    TitleHead   = ['Lilliput list: ']; stPlot.TitleHead = TitleHead;
    YLabel      = '% Correct';                  stPlot.YLabel = YLabel;
    XLabel      = '';                           stPlot.XLabel = XLabel;
    stPlot.ReverseData = 1; % from worst to best score
    stPlot.xLim = [0 2];
    stPlot.separation_series = 0.1;
    PlotMeans(stPlot, aceData, f0mData);
end

h = gcf;
set(gcf, 'PaperPositionMode','auto')
                
 
% idx_ACEown = find(data(:,column_strategy)==10); % baseline
% idx_ACE = find(data(:,column_strategy)==0);
% idx_F0m = find(data(:,column_strategy)==1);
% 
% aceDataOwn = data(idx_ACEown,[column_score_sentence column_score_word])*100;
% aceData = data(idx_ACE,[column_score_sentence column_score_word])*100; % in percentage. out of: 10 words
% f0mData = data(idx_F0m,[column_score_sentence column_score_word])*100;
% 
% column_scores_sentence_new  = 1; % after deleting
% column_scores_word_new      = 2; % after deleting
% 
% 
% aceDataOwn(:,column_scores_sentence_new)           = aceDataOwn(:,column_scores_sentence_new)/10;
% aceDataOwn(:,column_scores_word_new) = aceDataOwn(:,column_scores_word_new);
% aceData(:,column_scores_sentence_new)              = aceData(:,column_scores_sentence_new)/10;
% aceData(:,column_scores_word_new)  = aceData(:,column_scores_word_new);
% f0mData(:,column_scores_sentence_new)              = f0mData(:,column_scores_sentence_new)/10;
% f0mData(:,column_scores_word_new)  = f0mData(:,column_scores_word_new);
% 
% % Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% p.nStrategies       = 2;
% p.nSemitones        = 1;
% p.nReferences       = 2;
% nStrategies         = p.nStrategies;
% nSemitones          = p.nSemitones;
% nReferences         = p.nReferences;        stPlot.nPlots       = nReferences;
% 
% stPlot.nCols        = 2;
% stPlot.nRows        = 1;
% stPlot.xGap         = 15;
% stPlot.yGap         = 80;
% stPlot.fntsz        = 16;
% stPlot.markerSize   = 12;
% stPlot.margins      = [70 70 50 90];
% stPlot.figPos       = [0 0 1024 450];   % 225 each row
%                                         % left, right, top, bottom
% 
% stPlot.xTick = 1:nSemitones;
% stPlot.yTick = 0:20:100;
% stPlot.xTickLabel = {''};
% stPlot.YGrid = 'on';
% 
% SeriesLabel = {'ACE', 'F0mod'};
% stPlot.TitleSuffix = {' scores per word', ' scores per phoneme'};
% stPlot.SeriesLabel = SeriesLabel;
% TitleHead   = ['Lilliput list: ']; stPlot.TitleHead = TitleHead;
% YLabel      = '% Correct';                  stPlot.YLabel = YLabel;
% XLabel      = '';                           stPlot.XLabel = XLabel;
% stPlot.ReverseData = 1; % from worst to best score
% stPlot.xLim = [0 2];
% stPlot.separation_series = 0.1;
% 
% PlotMeans(stPlot, aceData, f0mData);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = gcf;
% set(gcf, 'PaperPositionMode','auto')

p.aceData = data(idx_ACE, [column_score_sentence column_score_word column_SNR]);
p.f0mData = data(idx_F0m, [column_score_sentence column_score_word column_SNR]);

p.column_score_sentence = 1;
p.column_score_word     = 2;
p.column_SNR            = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end