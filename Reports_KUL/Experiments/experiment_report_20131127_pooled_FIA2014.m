function experiment_report_20131127_pooled_FIA2014(options)
% function experiment_report_20131127_pooled_FIA2014(options)
%
%   Process CI Speech data
% 
% Previous data reports: 
%   2013 11 05: Holy Hour
%   2013 11 12: pooled data, release 2
%
% New data added:
%   MB, S12 on 2013 11 25, VlMatrix
%
% Pooled data session
%
% txt file with results:
%                                           nSubjects   SVN version
%       CI-Pooled-Sentence-LT-adaptive.txt  4 (11-14)   725
%                                           6 (11-16)
%       CI-Pooled-Word-LL.txt           
%
% Dependencies:
%   figLTScores_xPC
%   figVUScores_xPC
%   figSentenceScores_xPC
%   figLLScores_xPC
%   PlotMeans
%
% 3. Example:
%       opts.bDoLT = 1;
%       experiment_report_20131127_pooled_FIA2014(opts);
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2014
% Created on    : $$dd$$/$$mm$$/2014
% Last update on: 10/01/2014 (for Windows-TU/e compatibility) % Update this date manually
% Last use on   : 25/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
    bIsPaper = 0; % Calculations for Report
else
    bIsPaper = 1; % Calculations for IJA paper
end

session_date = '20131127';
subject_folder = 'ci-Pooled-report-3/'; % define folder where plots will be generated
options = Ensure_field(options,'bDoSpeechTests', 1);
options = Ensure_field(options,'bDoLT',1);
options = Ensure_field(options,'bDoStats',0);

bDoSpeechTests = options.bDoSpeechTests;
bDoLT       = options.bDoLT;
bDoStats    = options.bDoStats;

try
    info = getpaths;
    % info.TB_SP15
    % info.experiment_report  = [misc.svn.Meas 'Experiments/Experiment_reports/'];
    % info.result_folder      = [misc.svn.Meas 'Experiments/Results_XML/'];
    addpath(info.TB_SP15); % Ensure_field
catch
    warning('getpaths.m failed, probably file not found');
    % misc.svn.Meas = ['D:\SVN-KU-Leuven\alejandro\Meas' delim];
    misc.svn.Meas = [Get_TUe_data_paths('SVN_KUL') 'Meas' delim];
    info.experiment_report  = [misc.svn.Meas 'Experiments' delim 'Experiment_reports' delim];
    info.result_folder      = [misc.svn.Meas 'Experiments' delim 'Results_XML' delim];
end

% hoofd_folder    = [info.svn.Meas 'Experiments/'];

options.result_folder   = info.result_folder;
options     = Ensure_field(options,'dest_folder', [info.experiment_report session_date '-' subject_folder 'Figures/']);

% SubjectName = 'CIs';

file_LT_bef = [options.result_folder 'CI-Pooled-Sentence-LT-adaptive.txt']; % before training
file_LT     = [options.result_folder 'CI-Pooled-LT-after-training.txt'];
file_MT     = [options.result_folder 'CI-Pooled-MT-fixedSNR-after-training.txt']
file_LL     = [options.result_folder 'CI-Pooled-Word-LL.txt'];
 
options     = Ensure_field(options,'bSave',0);

if bDoSpeechTests == 1

    stPlot = [];
    
    if bDoLT
        
        subjects = [13 15 16 17]; % Everyone participated
        
        [xx pLT_bef]        = figLTScores_xPC(file_LT_bef,1);
        [xx pLT]            = figLTScores_xPC(file_LT    ,1);

        % Plotting LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pLT_bef = fullfil_subject_data(pLT_bef, subjects);
        pLT     = fullfil_subject_data(pLT    , subjects);
        
        clear stPlot
        
        stPlot.figPos       = [0 0 1024 300];
%         stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
%         stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
        stPlot.xTickLabel   = {'PM','DW','JBD','JS'};
        stPlot.xTickLabelPaper = {'S13','S15','S16','S17'};
        stPlot.xTick        = 1:length(stPlot.xTickLabel);
        stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
        stPlot.yLim         = [-8 8];
        stPlot.yTick        = [stPlot.yLim(1)+1:1:stPlot.yLim(2)-1];
        stPlot.SeriesLabel  = {'ACE','F0mod'};
        stPlot.SeriesColor  = {'w'  , 'k'};
%         stPlot.SeriesColor  = {[0.8 0.8 0.8], 'w'  , 'k'};
        stPlot.XLabel       = 'Subject';
        stPlot.YLabel       = 'SRT (dB)';
        stPlot.Title        = 'LIST-f sentences';
        
        ace2plot_bef    = [];
        f0m2plot_bef    = [];
        aceOwn2plot_bef = [];
        ace2plot        = [];
        f0m2plot        = [];
        aceOwn2plot     = [];
        
        for i = 1:length(subjects)
            idxSSA      = find(pLT_bef.aceData(:,pLT_bef.column_subject)==subjects(i) );
            ace2plot_bef= [ace2plot_bef pLT_bef.aceData(idxSSA,pLT_bef.column_score)];
            idxSSA      = find(   pLT.aceData(:,pLT.column_subject)==subjects(i) );
            ace2plot    = [ace2plot pLT.aceData(idxSSA,pLT.column_score)];

            idxSSB      = find(pLT_bef.f0mData(:,pLT_bef.column_subject)==subjects(i) );
            f0m2plot_bef= [f0m2plot_bef pLT_bef.f0mData(idxSSB,pLT_bef.column_score)];
            idxSSB      = find(pLT.f0mData(:,pLT.column_subject)==subjects(i) );
            f0m2plot    = [f0m2plot pLT.f0mData(idxSSB,pLT.column_score)];

            idxSSC      = find(pLT_bef.aceDataOwn(:,pLT_bef.column_subject)==subjects(i) );
            aceOwn2plot_bef = [aceOwn2plot_bef pLT_bef.aceDataOwn(idxSSC,pLT_bef.column_score)];
            idxSSC      = find(pLT.aceDataOwn(:,pLT.column_subject)==subjects(i) );
            aceOwn2plot = [aceOwn2plot pLT.aceDataOwn(idxSSC,pLT.column_score)];
        end
        
        if ~bIsPaper
            PlotMeans3series(stPlot, f0m2plot, ace2plot, aceOwn2plot);

            if bSave
                saveas(gcf,[dest_folder,'LIST-per-subject'],'epsc');
            end
            
            stPlot.bPlotIndividual = 0;
            stPlot.xTick        = [1];
            stPlot.xTickLabel   = {' '};

            stPlot.TitleHead    = [' '];
            stPlot.xLim         = [ 0   2];
            stPlot.XLabel       = '  ';
            stPlot.xTickLabel   = {'Pooled'};

            PlotMeans3series(stPlot,pLT.f0mData(:,2), pLT.aceData(:,2), pLT.aceDataOwn(:,2));  
            stPlot.TitleHead    = ['']; 
            xlim([ -2   2])
            hA = gca;
            set(hA,'YGrid','on')

            if bSave
                saveas(gcf,[dest_folder,'LIST-pooled'],'epsc');
            end
        end
                
        [mean2plot_bef std2plot_bef]= prepare_barplot(ace2plot_bef, f0m2plot_bef);
        [mean2plot     std2plot]    = prepare_barplot(ace2plot    , f0m2plot    );
        col_F0m_data        = 3;
        col_ACE_data        = 2;
        col_ACEown_data     = 1;
        
        % Colors = [0.70 0.70 0.70; 1 1 1;0 0 0];
        Colors = [0.7 0.7 0.7;1 1 1];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        stPlot.Title = 'LIST-f sentences';
        barweb(mean2plot_bef,std2plot_bef,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        grid on;
        tmp_h1 = gcf;
        set(tmp_h1,'Position', stPlot.figPos);

        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
        set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
        set(Handle,'YLim',stPlot.yLim)
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])

        set(tmp_h1, 'PaperPositionMode','auto')
        saveas(gcf,[Get_TUe_paths('outputs') 'LIST-f-4-subjects.eps']);

    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m s] = prepare_barplot(x, y, z)
% function [m s] = prepare_barplot(x, y, z)
%
% m - mean
% s - std

[mf0m       sf0m]       = Get_mean(x);
[mace       sace]       = Get_mean(y);
try
    [mace_own   sace_own]   = Get_mean(z);
end

[mf0m_tot,sf0m_tot]     = Get_mean(mf0m');
[mACE_tot,sACE_tot]     = Get_mean(mace');
try
    [mACEown_tot,sACEown_tot] = Get_mean(mace_own');
end
try
    m   = [mf0m' mace' mace_own'; nan(1,3); mf0m_tot mACE_tot mACEown_tot]; % NaN = trick to put data in a more beautiful way
    s   = [sf0m' sace' sace_own'; nan(1,3); sf0m_tot sACE_tot sACEown_tot];
catch
    m   = [mf0m' mace'          ; nan(1,2); mf0m_tot mACE_tot            ];
    s   = [sf0m' sace'          ; nan(1,2); sf0m_tot sACE_tot            ];
end
