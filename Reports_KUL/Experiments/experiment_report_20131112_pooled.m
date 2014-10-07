function experiment_report_20131112_pooled
% function experiment_report_20131112_pooled
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
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root_root       = '/home/tom/temp/alejandro/Meas/Experiments/'; 
dest_folder     = [root_root 'Experiment_reports/20131112-ci-Pooled-report-2/Figures/'];
hoofd_folder    = [root_root 'Results_XML/'];

SubjectName = 'CIs';

file_LIST_adaptive  = [hoofd_folder 'CI-Pooled-Sentence-LT-adaptive.txt'];
file_MT_fixedSNR    = [hoofd_folder 'CI-Pooled-MT-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'CI-Pooled-Word-LL.txt'];
 
bPlot               = 1;
bSave               = 1;
bDoSpeechTests      = 0;
bDoLB               = 0;
bDoPR               = 1;

if bDoSpeechTests == 1

    stPlot.bPlotIndividual  = 0;
    stPlot.YGrid            = 'on';
    stPlot.yLim             = [-7 7];
    stPlot.xLim             = [0 3];
    stPlot.xTickLabel       = {'Pooled'};
    [xx pLT]                = figLTScores_xPC(file_LIST_adaptive, bPlot, stPlot);
    xlim([-2 2])

    if bSave
        saveas(gcf,[dest_folder,'LT-adaptive-pooled'],'epsc');
    end
    
    try
        convertSpeech2SPSS_speech(pLT);
    end
    
    [h p]       = figLLScores_xPC(file_LL_fixedSNR, stPlot);
    close;
    
    try
        convertLL2SPSS_speech(p);
    end
    
    % Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idx_QuietA = find(p.aceData(:,p.column_SNR)==99);
    idx_QuietB = find(p.f0mData(:,p.column_SNR)==99);
    idx_QuietC = find(p.aceDataOwn(:,p.column_SNR)==99);
    idx_NoiseA = find(p.aceData(:,p.column_SNR)==10);
    idx_NoiseB = find(p.f0mData(:,p.column_SNR)==10);
    idx_NoiseC = find(p.aceDataOwn(:,p.column_SNR)==10);

    stPlot.nCols        = 1;
    stPlot.nRows        = 1;
    stPlot.xGap         = 15;
    stPlot.yGap         = 80;
    stPlot.fntsz        = 16;
    stPlot.markerSize   = 12;
    stPlot.margins      = [70 70 50 90];
    stPlot.figPos       = [0 0 1024 450]; % 225 each row % left, right, top, bottom
    stPlot.xLim         = [  0  3];
    stPlot.yLim         = [15 105];
    stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
    stPlot.xTickLabel   = {'Q','+10'};
    stPlot.YGrid        = 'on';
    stPlot.TitleSuffix  = {''};
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE','ACE baseline'};
    stPlot.XLabel       = 'SNR (dB)';
    stPlot.ReverseData  = 0; % from worst to best score
    stPlot.separation_series = 0.1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lilliput Scores: 

    stPlot.figPos       = [0 0 1024 450];
    stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
    stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
    stPlot.xTick        = 1:length(stPlot.xTickLabel);
    stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
    stPlot.XLabel       = 'Subject';
    stPlot.YLabel       = 'Scores (%)';
    MaxScore_p          = 30;
    stPlot.bPlotIndividual = 1;

    subjects = [11 12 13 14 15 16 17];

    SNR = [99 10];
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
            stPlot.TitleHead    = ['Lilliput in quiet (Results for ' SubjectName ')'];
        else
            stPlot.TitleHead    = ['Lilliput, SNR = 10 dB (Results for ' SubjectName ')']; 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [mf0m sf0m]             = Get_mean(f0m2plot);
        [mace sace]             = Get_mean(ace2plot);
        [mace_own sace_own]     = Get_mean(aceOwn2plot);
        
        [mf0m_tot,sf0m_tot]     = Get_mean(mf0m');
        [mACE_tot,sACE_tot]     = Get_mean(mace');
        [mACEown_tot,sACEown_tot] = Get_mean(mace_own');
        
        mean2plot   = [mf0m' mace' mace_own'; nan(1,3); mf0m_tot mACE_tot mACEown_tot];
        std2plot    = [sf0m' sace' sace_own'; nan(1,3); sf0m_tot sACE_tot sACEown_tot];
        
        Colors = [0 0 0; 1 1 1; 0.70 0.70 0.70];
        figure
        barweb(mean2plot,std2plot,[],[],stPlot.TitleHead,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        tmp_h2 = gcf;
        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
        set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
        set(Handle,'YLim',[10 110])
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if SNR(k) == 99
            if bSave
                %saveas(tmp_h1,[dest_folder,'Lilliput-p-scores-quiet-per-subject'],'epsc');
                saveas(tmp_h2,[dest_folder,'Lilliput-p-scores-quiet-per-subject-bar'],'epsc');
            end
        end

        if SNR(k) == 10
            if bSave
                %saveas(tmp_h1,[dest_folder,'Lilliput-p-scores-noise-per-subject'],'epsc');
                saveas(tmp_h2,[dest_folder,'Lilliput-p-scores-noise-per-subject-bar'],'epsc');
            end
        end

    end
    
    for i = 1:length(subjects) % Individual pooled results. One plot per subject
        stHandle = quick_LL_one_subject(p, subjects(i));
 
        if bSave
            saveas(stHandle.h1     ,[dest_folder,'Lilliput-one-subject-report-S'   ,num2str(subjects(i))],'epsc');
            close
        end
    end

    stPlot.bPlotIndividual = 0;
    stPlot.xTick        = [1];
    stPlot.xTickLabel   = {' '};

    stPlot.TitleHead    = [' '];
    stPlot.xLim         = [ 0   2];
    stPlot.xTickLabel   = 'Pooled';
    stPlot.XLabel       = ' ';

%     PlotMeans3series(stPlot,[p.f0mData(idx_QuietB,p.column_score_phoneme   )/MaxScore_p]*100, ...
%                             [p.aceData(idx_QuietA ,p.column_score_phoneme  )/MaxScore_p]*100, ...
%                             [p.aceDataOwn(idx_QuietC,p.column_score_phoneme)/MaxScore_p]*100);  
    xlim([-2 2])

    stPlot.TitleHead    = ['']; 
    stPlot.xTickLabel   = {' '};

%     if bSave
%         saveas(gcf,[dest_folder,'Lilliput-p-scores-quiet-pooled'],'epsc');
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stPlot.TitleHead    = ['Lilliput, SNR = 10 dB (Results for ' SubjectName ')']; 
    stPlot.line_color   = 'w'; % no line
    stPlot.xLim         = [ 0   3+2.5];
    stPlot.XLabel       = 'Subject';
    stPlot.YLabel       = 'Scores (%)';
    MaxScore_v          = 10;
    MaxScore_c          = 20;
    MaxScore_p          = 30;
    MaxScore_w          = 10;

    stPlot.bPlotIndividual = 0;
    stPlot.xTick        = [1];
    stPlot.xTickLabel   = {' '};

    stPlot.TitleHead    = [' '];
    stPlot.xLim         = [0   2];
    stPlot.xTickLabel   = 'Pooled';
    stPlot.XLabel       = ' ';

%     PlotMeans3series(stPlot,[p.f0mData(idx_NoiseB,p.column_score_phoneme   )/MaxScore_p]*100, ...
%                             [p.aceData(idx_NoiseA ,p.column_score_phoneme  )/MaxScore_p]*100, ...
%                             [p.aceDataOwn(idx_NoiseC,p.column_score_phoneme)/MaxScore_p]*100);  
%     xlim([-2 2])

    stPlot.TitleHead    = ['']; 

%     if bSave
%         saveas(gcf,[dest_folder,'Lilliput-p-scores-noise-pooled'],'epsc');
%     end

    % Plotting LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pLT = fullfil_subject_data(pLT, subjects);

    % stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
    stPlot.xTickLabel   = {'S11','S12','S13','S14','S15','S16','S17'};
    stPlot.xTick        = 1:length(stPlot.xTickLabel);
    stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
    stPlot.XLabel       = 'Subject';
    stPlot.YLabel       = 'SRT (dB)';

    stPlot.yTick            = [-7+1:2:7];

    stPlot.bPlotIndividual  = 1;
    stPlot.yLim             = [-7   7];

    pLT.column_SRT = 2;

    ace2plot        = [];
    f0m2plot        = [];
    aceOwn2plot     = [];
    for i = 1:length(subjects)
        idxSSA      = find(   pLT.aceData(:,pLT.column_subject)==subjects(i) );
        ace2plot    = [ace2plot pLT.aceData(idxSSA,pLT.column_SRT)];

        idxSSB      = find(   pLT.f0mData(:,pLT.column_subject)==subjects(i) );
        f0m2plot    = [f0m2plot pLT.f0mData(idxSSB,pLT.column_SRT)];

        idxSSC      = find(pLT.aceDataOwn(:,pLT.column_subject)==subjects(i) );
        aceOwn2plot = [aceOwn2plot pLT.aceDataOwn(idxSSC,pLT.column_SRT)];
    end

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

    if bSave
        saveas(gcf,[dest_folder,'LIST-pooled'],'epsc');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATRIX: 

    stPlot.bPlotIndividual = 1;
    [handleFig pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);

    if bSave
        saveas(handleFig(1),[dest_folder,'MT-SNR10dB-pooled'],'epsc');
        saveas(handleFig(2),[dest_folder,'MT-SNR05dB-pooled'],'epsc');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix word-scores: 

    stPlot.figPos       = [0 0 1024 450];
    stPlot.xTickLabel   = {'RP','MB','PM','JL','DW','JBD','JS'};
    %stPlot.xTickLabelPaper = {'S11','S12','S13','S14','S15','S16','S17'};
    stPlot.xTick        = 1:length(stPlot.xTickLabel);
    stPlot.xLim         = [ 0   max(stPlot.xTick)+2];
    stPlot.XLabel       = 'Subject';
    stPlot.YLabel       = 'Word scores (%)';
    stPlot.bPlotIndividual = 1;
    stPlot.yLim         = [15 105];
    stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
    stPlot.YGrid        = 'on';
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};
    stPlot.ReverseData  = 1; % from worst to best score
    stPlot.separation_series = 0.1;

    SNR = [10 5];
    
    convertMT2SPSS_speech(pMT);
    
    for k = 1:2

        ace2plot        = [];
        f0m2plot        = [];
        for i = 1:length(subjects)
            idxSSA      = find(   pMT.aceData(:,pMT.column_subject)==subjects(i) &    pMT.aceData(:,pMT.column_SNR)==SNR(k));
            ace2plot    = [ace2plot pMT.aceData(idxSSA,pMT.column_score_word)];

            idxSSB      = find(   pMT.f0mData(:,pMT.column_subject)==subjects(i) &    pMT.f0mData(:,pMT.column_SNR)==SNR(k));
            f0m2plot    = [f0m2plot pMT.f0mData(idxSSB,pMT.column_score_word)];
        end

        PlotMeans(stPlot, f0m2plot, ace2plot);
        tmp_h1 = gcf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [mf0m sf0m]             = Get_mean(f0m2plot);
        [mace sace]             = Get_mean(ace2plot);
        [mf0m_tot,sf0m_tot]     = Get_mean(mf0m');
        [mACE_tot,sACE_tot]     = Get_mean(mace');

        mean2plot   = [mf0m' mace'; nan(1,2); mf0m_tot mACE_tot];
        std2plot    = [sf0m' sace'; nan(1,2); sf0m_tot sACE_tot];
        
        Colors = [0 0 0; 1 1 1];
        
        figure
        barweb(mean2plot,std2plot,[],[],stPlot.TitleHead,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        tmp_h2 = gcf;
        Handle = gca;
        set(Handle,'YTick',stPlot.yTick)
        set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
        set(Handle,'XLim',[0 max(stPlot.xTick)+2+0.5])
        set(Handle,'YLim',[10 110])
        set(Handle,'XTickLabel',[stPlot.xTickLabelPaper,'AVG'])

        if SNR(k) == 10
            if bSave
                saveas(tmp_h1,[dest_folder,'Vlaamse-Matrix-w-scores-SNRp10'],'epsc');
                saveas(tmp_h2,[dest_folder,'Vlaamse-Matrix-w-scores-SNRp10-bar'],'epsc');
            end
        end

        if SNR(k) == 5
            if bSave
                saveas(tmp_h1,[dest_folder,'Vlaamse-Matrix-w-scores-SNRp05'],'epsc');
                saveas(tmp_h2,[dest_folder,'Vlaamse-Matrix-w-scores-SNRp05-bar'],'epsc');
            end
        end

    end
end

if bDoPR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loudness Balancing

    stPlot = [];
    stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};

    if bDoLB
        % When required type: 'ci-'
        [deltaACE, deltaF0mod]  = figLBScores_xPC(hoofd_folder,{'UW'}, bPlot);

        if bSave
            Min4plot = round(min(min([deltaACE; deltaF0mod]))-5);
            Max4plot = round(max(max([deltaACE; deltaF0mod]))+5);
            ylim([Min4plot Max4plot])

            saveas(gcf,[dest_folder,'LB-results'],'epsc');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pitch ranking

    stSubject.firstname = 'ci-';
    stSubject.truncate_data = 0;

    [aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
    % To confirm the data, type 1
    % Exclude Maria Brughmans' data

    convertPR2SPSS(aceData,f0mData);

    if bSave
        saveas(gcf,[dest_folder,'PR-results-pooled'],'epsc');
    end

    stSubject.firstname = 'ci-Romain';
    % each time new tests are stored, select these values
    [aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-RP'],'epsc');
    end

    stSubject.firstname = 'ci-Patrick';
    [aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-PM'],'epsc');
    end

    stSubject.firstname = 'ci-Jan';
    [aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-JL'],'epsc');
    end

    stSubject.firstname = 'ci-Jean-Baptist';
    [aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
    if bSave
        saveas(gcf,[dest_folder,'PR-results-JBD'],'epsc');
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end