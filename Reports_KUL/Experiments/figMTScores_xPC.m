function [h p] = figMTScores_xPC(filename,stPlot,nSentences)
% function [h p] = figMTScores_xPC(filename,stPlot,nSentences)
%
% filename - txt file containing VlMatrix results
%
% column_report:    1 = reported before 31-10-2013
% Strategy = 0 = ACE
% Strategy = 1 = F0mod
%
% % Example:
%        stPlot.bPlotIndividual = 1;
%        main_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/'; 
%        file_MT = [main_folder 'ci-Maria_Brughmans/20131125-MT/20131125-CI-MB-VlMatrix.txt';
%        [xxx pMT] = figMTScores_xPC(file_MT,stPlot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    nSentences = 20; % doublelist
end
if nargin < 2
    bPlot = 0;
else
    bPlot = 1;
end

% List_of_files = {   'confidenceInterval.m', ...
%                     'get_significance_level.m', ...
%                     'setupPlotConf.m'};
%                 
% [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

column_subject          = 1;
column_SNR              = 2;
column_strategy         = 3; 
column_score_sentence   = column_strategy + 1;
column_score_name       = column_strategy + 3;
column_score_verb       = column_strategy + 4;
column_score_numeral    = column_strategy + 5;
column_score_color      = column_strategy + 6;
column_score_object     = column_strategy + 7;
column_num_list         = column_strategy + 2;   
column_nSentences       = column_strategy + 8;
column_report           = column_strategy + 9;

import_scores_from_txt(filename);
data = evalin('base','data');
textdata = evalin('base','textdata');

idx_ACEown = find(data(:,column_strategy)==10); % baseline
idx_ACE = find(data(:,column_strategy)==0);
idx_F0m = find(data(:,column_strategy)==1);

try
    aceDataOwn = data(idx_ACEown,[column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object])./data(idx_ACEown,column_nSentences)*100;
catch
    aceDataOwn = data(idx_ACEown,[column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object])*100;
end
aceData = data(idx_ACE      ,[column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object])./repmat(data(idx_ACE,column_nSentences),1,6)*100; % in percentage. out of: 10 words
f0mData = data(idx_F0m      ,[column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object])./repmat(data(idx_F0m,column_nSentences),1,6)*100;

aceData_SNR = data(idx_ACE,2);
f0mData_SNR = data(idx_F0m,2);
% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.nStrategies       = 2;
p.nSemitones        = 1;

stPlot.nCols        = 1;
stPlot.nRows        = 1;
stPlot.xGap         = 15;
stPlot.yGap         = 80;
stPlot.fntsz        = 16;
stPlot.markerSize   = 12;
stPlot.margins      = [70 70 50 90];
stPlot.figPos       = [0 0 1024 450];   % 225 each row
                                        % left, right, top, bottom

stPlot.xTick    = 1:size(aceData,2);
stPlot.yTick    = 0:20:100;
stPlot.xTickLabel = {''};
stPlot.YGrid    = 'on';

stPlot.SeriesLabel = {'xPC F0mod', 'xPC ACE'};
stPlot.TitleHead = 'Flemish Matrix, SNR = 10 dB';
stPlot.YLabel   = '% Correct';
stPlot.XLabel   = {'Score type'};
stPlot.ReverseData = 1; % from worst to best score
stPlot.xLim     = [0    8];
stPlot.yLim     = [0 120];
stPlot.xTickLabel = {'sentence', 'name','verb','numeral','color','object'};
stPlot.separation_series = 0.1;
stPlot.line_color = 'w';
stPlot.ReverseData = 0;

idx_SNRA = find(aceData_SNR==10);
idx_SNRB = find(aceData_SNR==10);

if bPlot
    PlotMeans(stPlot, f0mData(idx_SNRB,:), aceData(idx_SNRA,:));
    h(end+1) = gcf;
    if strcmp(stPlot.YGrid, 'on')
        set(gca,'YGrid','on')
    end
    set(gcf, 'PaperPositionMode','auto')

    stPlot.TitleHead = ['Flemish Matrix, SNR = 5 dB'];
    idx_SNRA = find(aceData_SNR==5);
    idx_SNRB = find(aceData_SNR==5);
    PlotMeans(stPlot, f0mData(idx_SNRB,:), aceData(idx_SNRA,:));
    h(end+1) = gcf;
    if strcmp(stPlot.YGrid, 'on')
        set(gca,'YGrid','on')
    end
    set(gcf, 'PaperPositionMode','auto')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meanACE = (sum( data(idx_ACE, [column_score_name column_score_verb column_score_numeral column_score_color column_score_object])' )./data(idx_ACE,column_nSentences)' )/5;
meanACE = meanACE'*100;

meanF0m = (sum( data(idx_F0m, [column_score_name column_score_verb column_score_numeral column_score_color column_score_object])' )./data(idx_F0m,column_nSentences)' )/5;
meanF0m = meanF0m'*100;


if ~exist('column_SNR','var')
    p.aceData = [data(idx_ACE, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject]) meanACE data(idx_ACE,column_report)];
    p.f0mData = [data(idx_F0m, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject]) meanF0m data(idx_F0m,column_report)];
    
    p.column_score_word = 8;
    
    if length(idx_ACEown) ~= 0
        p.aceDataOwn = data(idx_ACEown, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject]);
    end
else
    p.aceData = [data(idx_ACE, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject column_SNR]) meanACE data(idx_ACE,column_report)];
    p.f0mData = [data(idx_F0m, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject column_SNR]) meanF0m data(idx_F0m,column_report)];
    
    p.column_SNR        = 8;
    p.column_score_word = 9;
    
    if length(idx_ACEown) ~= 0
        p.aceDataOwn = data(idx_ACEown, [column_score_sentence column_score_name column_score_verb column_score_numeral column_score_color column_score_object column_subject column_SNR column_report]);
    end
    
end
p.column_report         = p.column_score_word+1;
p.column_score_sentence = 1;
p.column_score_name     = 2;
p.column_score_verb     = 3;
p.column_score_numeral  = 4;
p.column_score_color    = 5;
p.column_score_object   = 6;
p.column_subject        = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
