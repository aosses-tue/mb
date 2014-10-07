function experiment_report_20131127_pooled(options)
% function experiment_report_20131127_pooled(options)
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
options = Ensure_field(options,'bDoLL',1);
options = Ensure_field(options,'bDoMT',1);
options = Ensure_field(options,'bDoStats',0);

bDoSpeechTests = options.bDoSpeechTests;
bDoLT       = options.bDoLT;
bDoLL       = options.bDoLL;
bDoMT       = options.bDoMT;
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
        
        subjects = [11 12 13 14 15 16 17]; % Everyone participated
        
        [xx pLT_bef]        = figLTScores_xPC(file_LT_bef,1);
        [xx pLT]            = figLTScores_xPC(file_LT    ,1);
        
        xlim([-2 2])

        if options.bSave && ~bIsPaper
            saveas(gcf,[dest_folder,'LIST-pooled'],'epsc');
        end
    
        if bDoStats
            try
                convertSpeech2SPSS_speech(pLT_bef);
                convertSpeech2SPSS_speech(pLT    );
                % Manually copied to Meas/Meas/Experiments/Results_XML/Data4SPSS/data4SPSS-LT.txt
            end
        end

        % Plotting LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pLT_bef = fullfil_subject_data(pLT_bef, subjects);
        pLT     = fullfil_subject_data(pLT    , subjects);
        
        clear stPlot
        
        stPlot.figPos       = [0 0 1024 300];
        stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
        stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
        stPlot.xTick        = 1:length(stPlot.xTickLabel);
        stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
        stPlot.yLim         = [-8 8];
        stPlot.yTick        = [stPlot.yLim(1)+1:1:stPlot.yLim(2)-1];
        stPlot.SeriesLabel  = {'ACE OP', 'ACE xPC','F0mod'};
        stPlot.SeriesColor  = {[0.8 0.8 0.8], 'w'  , 'k'};
        stPlot.XLabel       = 'Subject';
        stPlot.YLabel       = 'SRT (dB)';
        stPlot.Title        = 'LIST sentences';
        
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
                
        [mean2plot_bef std2plot_bef]= prepare_barplot(aceOwn2plot_bef, ace2plot_bef, f0m2plot_bef);
        [mean2plot     std2plot]    = prepare_barplot(aceOwn2plot    , ace2plot    , f0m2plot    );
        col_F0m_data        = 3;
        col_ACE_data        = 2;
        col_ACEown_data     = 1;
        
        Colors = [0.70 0.70 0.70; 1 1 1;0 0 0];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        stPlot.Title = 'LIST sentences';
        barweb(mean2plot_bef,std2plot_bef,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        
        tmp_h1 = gcf;
        set(tmp_h1,'Position', stPlot.figPos);

        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
        set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
        set(Handle,'YLim',stPlot.yLim)
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])

        set(tmp_h1, 'PaperPositionMode','auto')
        
        %------------------------------------------------------------------
        figure
        stPlot.Title = 'LIST sentences, after training';
        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        
        tmp_h2 = gcf;
        set(tmp_h2,'Position', stPlot.figPos);

        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
        set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
        set(Handle,'YLim',stPlot.yLim)
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])

        set(tmp_h2, 'PaperPositionMode','auto')
        
        Colors = [1 1 1; 0.35 0.35 0.35; 0 0 0];
        
        mTraining = [mean2plot(2, col_ACE_data) mean2plot_bef(2, col_F0m_data) mean2plot(2, col_F0m_data); ...
                     mean2plot(4, col_ACE_data) mean2plot_bef(4, col_F0m_data) mean2plot(4, col_F0m_data); ...
                     mean2plot(7, col_ACE_data) mean2plot_bef(7, col_F0m_data) mean2plot(7, col_F0m_data)];
        sTraining = [ std2plot(2, col_ACE_data)  std2plot_bef(2, col_F0m_data)  std2plot(2, col_F0m_data) ; ...
                      std2plot(4, col_ACE_data)  std2plot_bef(4, col_F0m_data)  std2plot(4, col_F0m_data) ; ...
                      std2plot(7, col_ACE_data)  std2plot_bef(7, col_F0m_data)  std2plot(7, col_F0m_data)];
        
        stPlot.figPos       = [0 0 1024 300];
        stPlot.SeriesLabel  = {'ACE xPC', 'F0mod (before training)','F0mod (after training)'};
        stPlot.xTickLabelPaper = {'S12', 'S14', 'S17'};
        stPlot.xTickLabel   = {'MB','JL','JS'};
        stPlot.xTick        = 1:length(stPlot.xTickLabel);
        stPlot.xLim         = [ 0.5   max(stPlot.xTick)+0.5];
        stPlot.yLim         = [-4 9];
        stPlot.yTick        = [stPlot.yLim(1)+1:1:stPlot.yLim(2)-1];
        
        stPlot.XLabel       = 'Subject';
        stPlot.YLabel       = 'SRT (dB)';
        figure
        stPlot.Title = 'LIST sentences for subjects S12, S14 and S17';
        barweb(mTraining,sTraining,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        
        tmp_h3 = gcf;
        set(tmp_h3,'Position', stPlot.figPos);

        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',stPlot.xTick)
        set(Handle,'XLim',[0 max(stPlot.xTick)+1])
        set(Handle,'YLim',stPlot.yLim)
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper])

        set(tmp_h3, 'PaperPositionMode','auto')
        
        % End of LIST-f figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        if options.bSave
            saveas(tmp_h1,[options.dest_folder,'list-before-training'],'epsc');
            saveas(tmp_h2,[options.dest_folder,'list-after-training'] ,'epsc');
            saveas(tmp_h3,[options.dest_folder,'list-bef-after-training'] ,'epsc');
            disp([mfilename '.m: Figures for LIST-f saved in ' options.dest_folder])
        end
        
    end
    % End of LIST-f processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if bDoLL
        
        subjects    = [11 12 13 14 15 16 17]; % Everyone participated
        SNR         = [99 10]; % Experiments conducted in quiet and in noise
        [h p]       = figLLScores_xPC(file_LL);
        close;
       
        try
            convertLL2SPSS_speech(p);
        end
    
        % Plotting: Lilliput Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear stPlot
        
        stPlot.figPos       = [0 0 1024 300];
        stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
        stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
        stPlot.xTick        = 1:length(stPlot.xTickLabel);
        stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
        stPlot.yLim         = [5 105];
        stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
        stPlot.SeriesLabel  = {'ACE OP', 'ACE xPC','F0mod'};
        stPlot.XLabel       = 'Subject';
        stPlot.YLabel       = 'Phoneme score (%)';
        MaxScore_p          = 30;

        for k = 1:2

            ace2plot        = [];
            f0m2plot        = [];
            aceOwn2plot     = [];
            for i = 1:length(subjects)
                idxSSA      = find(   p.aceData(:,p.column_subject)==subjects(i) &    p.aceData(:,p.column_SNR)==SNR(k));
                ace2plot    = [ace2plot p.aceData(idxSSA,p.column_score_phoneme)/MaxScore_p*100];

                idxSSB      = find(   p.f0mData(:,p.column_subject)==subjects(i) &    p.f0mData(:,p.column_SNR)==SNR(k));
                f0m2plot    = [f0m2plot p.f0mData(idxSSB,p.column_score_phoneme)/MaxScore_p*100];

                idxSSC      = find(p.aceDataOwn(:,p.column_subject)==subjects(i) & p.aceDataOwn(:,p.column_SNR)==SNR(k));
                aceOwn2plot = [aceOwn2plot p.aceDataOwn(idxSSC,p.column_score_phoneme)/MaxScore_p*100];
            end

            if SNR(k) == 99
                stPlot.Title    = ['Lilliput words in quiet'];
            else
                stPlot.Title    = ['Lilliput words in noise, SNR = ' num2str(SNR(k)) ' dB']; 
            end
            
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [mean2plot std2plot] = prepare_barplot(aceOwn2plot, ace2plot, f0m2plot);
            
            Colors = [0.7 0.7 0.7; 1 1 1; 0 0 0];
            
            figure
            barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);

            tmp_h2 = gcf;
            set(tmp_h2,'Position', stPlot.figPos);

            Handle = gca;
            set(Handle,'YTick',stPlot.yTick)
            set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
            set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
            set(Handle,'YLim',[5 110])
            set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])

            set(tmp_h2, 'PaperPositionMode','auto')
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if SNR(k) == 99
                if options.bSave
                    saveas(tmp_h2,[options.dest_folder,'lilliputQ-new'],'epsc');
                end
            end
    
            if SNR(k) == 10
                if options.bSave
                    saveas(tmp_h2,[options.dest_folder, 'lilliput10dB-new'],'epsc');
                    disp([mfilename '.m: Figures for Lilliput saved in ' options.dest_folder])
                end
            end
        end
        
    end
    % End of Lilliput processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if bDoMT % Flemish MATRIX:
        
        clear stPlot
        
        stPlot.bPlotIndividual = 1;
        [handleFig pMT] = figMTScores_xPC(file_MT);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Matrix word-scores: 

        stPlot.figPos       = [0 0 1024 300];
        stPlot.xTickLabel   = {'RP','MB','PM','JL','WD','JBD','JS'};
        stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
        stPlot.xTick        = 1:length(stPlot.xTickLabel);
        stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
        stPlot.XLabel       = 'Subject';
        stPlot.YLabel       = 'Word scores (%)';
        stPlot.yLim         = [5 105];
        stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
        stPlot.YGrid        = 'on';
        stPlot.SeriesLabel  = {'ACE', 'F0mod'};
        stPlot.separation_series = 0.1;

        SNR = [10 5];

        try
            convertMT2SPSS_speech(pMT);
        end

        for k = 1 % for k = 1:2 % We will do just SNR = 10 dB

            subjects = [11 12 13 14 15 16 17]; % Everyone participated
            
            ace2plot        = [];
            f0m2plot        = [];
            
            for i = 1:length(subjects)
                idxSSA      = find(   pMT.aceData(:,pMT.column_subject)==subjects(i) &    pMT.aceData(:,pMT.column_SNR)==SNR(k));
                ace2plot    = [ace2plot pMT.aceData(idxSSA,pMT.column_score_word)];
    
                idxSSB      = find(   pMT.f0mData(:,pMT.column_subject)==subjects(i) &    pMT.f0mData(:,pMT.column_SNR)==SNR(k));
                f0m2plot    = [f0m2plot pMT.f0mData(idxSSB,pMT.column_score_word)];
            end

            if SNR(k) == 10
                stPlot.Title    = ['Flemish Matrix sentences, SNR = 10 dB'];
            else
                stPlot.Title    = ['Flemish Matrix sentences, SNR = 5 dB']; 
            end
         
            [mean2plot std2plot] = prepare_barplot(ace2plot, f0m2plot);
            Colors = [1 1 1; 0 0 0];

            figure
            barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
            tmp_h2 = gcf;
            Handle = gca;
            set(Handle,'YTick',stPlot.yTick)
            set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
            set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
            set(Handle,'YLim',stPlot.yLim)
            set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])
            set(tmp_h2,'Position',stPlot.figPos)

            set(tmp_h2, 'PaperPositionMode','auto')
            
            if options.bSave
                figname = [options.dest_folder,'FlemishMatrix',num2str(SNR(k)),'dB'];
                saveas(tmp_h2,figname,'epsc');
                disp([mfilename '.m: Figures for VlMatrix saved in ' options.dest_folder])
                disp([mfilename '.m: Figure saved as ' figname]);
            end
           
        end

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
