function experiment_report_20131004
% function experiment_report_20131004
%
% Session with Jan Leys
%
% Dependencies:
%   figLTScores_xPC
%   figVUScores_xPC
%   figSentenceScores_xPC
%   figLLScores_xPC
%   PlotMeans
%
% Programmed by Alejandro, adapt to actual Patrick Meul's data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
dest_folder     = [hoofd_folder 'ci-pooled/Figures/'];

SubjectName = 'CI 1,3,4';

file_LIST_adaptive  = [hoofd_folder 'CI-Pooled-Sentence-LT-adaptive.txt'];
% file_VU_adaptive    = [hoofd_folder 'CI-Patrick-VU-adaptive.txt'];
% file_MT_fixedSNR    = [hoofd_folder 'CI-Patrick-MT-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'CI-Pooled-Word-LL.txt'];
 
bPlot = 1;
bSave = 1;

stPlot.bPlotIndividual  = 0;
stPlot.YGrid            = 'on';
stPlot.yLim             = [-7 7];
stPlot.xLim             = [0 3];

[xx pLT]                = figLTScores_xPC(file_LIST_adaptive, bPlot, stPlot);

if bSave
    saveas(gcf,[dest_folder,'LT-adaptive-pooled'],'epsc');
end

% stPlot.yLim = [0 10];
% [xx pVUadaptive] = figVUScores_xPC(file_VU_adaptive, bPlot, stPlot);
% if bSave
%     saveas(gcf,[dest_folder,'VU-adaptive-results'],'epsc');
% end

[h p]       = figLLScores_xPC(file_LL_fixedSNR, stPlot);
close;


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
stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod','ACE baseline'};
stPlot.TitleHead    = ['Lilliput material (Results for CI 3 - Patrick Meul)']; 
stPlot.XLabel       = 'SNR (dB)';
stPlot.ReverseData  = 0; % from worst to best score
stPlot.separation_series = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lilliput Scores: 

stPlot.TitleHead    = ['Lilliput in quiet (Results for ' SubjectName ')']; 
stPlot.figPos       = [0 0 1024 450];
stPlot.xLim         = [ 0   5.5];
stPlot.xTick        = [1 2 3];
stPlot.xTickLabel   = {'RP','PM', 'JL'};
stPlot.XLabel       = 'Subject';
stPlot.YLabel       = 'Scores (%)';
% score2plot          = p.column_score_vowel;
MaxScore_v          = 10;
MaxScore_c          = 20;
MaxScore_p          = 30;
MaxScore_w          = 10;
% idxS11A = find(p.aceData(idx_QuietA,p.column_subject)==11);
% idxS11B = find(p.f0mData(idx_QuietB,p.column_subject)==11);
% idxS11C = find(p.aceDataOwn(idx_QuietC,p.column_subject)==11);
% idxS13A = find(p.aceData(idx_QuietA,p.column_subject)==13);
% idxS13B = find(p.f0mData(idx_QuietB,p.column_subject)==13);
% idxS13C = find(p.aceDataOwn(idx_QuietC,p.column_subject)==13);
% idxS14A = find(p.aceData(idx_QuietA,p.column_subject)==14);
% idxS14B = find(p.f0mData(idx_QuietB,p.column_subject)==14);
% idxS14C = find(p.aceDataOwn(idx_QuietC,p.column_subject)==14);

idxS11A = find(   p.aceData(:,p.column_subject)==11 &    p.aceData(:,p.column_SNR)==99);
idxS11B = find(   p.f0mData(:,p.column_subject)==11 &    p.f0mData(:,p.column_SNR)==99);
idxS11C = find(p.aceDataOwn(:,p.column_subject)==11 & p.aceDataOwn(:,p.column_SNR)==99);
idxS13A = find(   p.aceData(:,p.column_subject)==13 &    p.aceData(:,p.column_SNR)==99);
idxS13B = find(   p.f0mData(:,p.column_subject)==13 &    p.f0mData(:,p.column_SNR)==99);
idxS13C = find(p.aceDataOwn(:,p.column_subject)==13 & p.aceDataOwn(:,p.column_SNR)==99);
idxS14A = find(   p.aceData(:,p.column_subject)==14 &    p.aceData(:,p.column_SNR)==99);
idxS14B = find(   p.f0mData(:,p.column_subject)==14 &    p.f0mData(:,p.column_SNR)==99);
idxS14C = find(p.aceDataOwn(:,p.column_subject)==14 & p.aceDataOwn(:,p.column_SNR)==99);



PlotMeans3series(stPlot,[p.aceData(idxS11A,p.column_score_phoneme   )/MaxScore_p  p.aceData(idxS13A,p.column_score_phoneme   )/MaxScore_p p.aceData(idxS14A,p.column_score_phoneme   )/MaxScore_p]*100, ...
                        [p.f0mData(idxS11B,p.column_score_phoneme   )/MaxScore_p  p.f0mData(idxS13B,p.column_score_phoneme   )/MaxScore_p p.f0mData(idxS14B,p.column_score_phoneme   )/MaxScore_p]*100, ...
                        [p.aceDataOwn(idxS11C,p.column_score_phoneme)/MaxScore_p  p.aceDataOwn(idxS13C,p.column_score_phoneme)/MaxScore_p p.aceDataOwn(idxS14C,p.column_score_phoneme)/MaxScore_p]*100);  
if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-quiet-per-subject'],'epsc');
end

stPlot.bPlotIndividual = 0;
stPlot.xTick        = [1];
stPlot.xTickLabel   = {' '};

stPlot.TitleHead    = [' '];
stPlot.xLim         = [ 0   2];
stPlot.XLabel       = 'Pooled';

PlotMeans3series(stPlot,[p.aceData(idx_QuietA ,p.column_score_phoneme   )/MaxScore_p]*100, ...
                        [p.f0mData(idx_QuietB,p.column_score_phoneme   )/MaxScore_p]*100, ...
                        [p.aceDataOwn(idx_QuietC,p.column_score_phoneme)/MaxScore_p]*100);  
stPlot.TitleHead    = ['']; 
stPlot.xTickLabel   = {' '};

if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-quiet-pooled'],'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stPlot.TitleHead    = ['Lilliput, SNR = 10 dB (Results for ' SubjectName ')']; 
stPlot.line_color   = 'w'; % no line
stPlot.bPlotIndividual = 1;
stPlot.xTick        = [1 2 3];
stPlot.xLim         = [ 0   3+2.5];
stPlot.xTickLabel   = {'RP','PM', 'JL'};
stPlot.XLabel       = 'Subject';
stPlot.YLabel       = 'Scores (%)';
MaxScore_v          = 10;
MaxScore_c          = 20;
MaxScore_p          = 30;
MaxScore_w          = 10;

idxS11A = find(   p.aceData(:,p.column_subject)==11 &    p.aceData(:,p.column_SNR)==10);
idxS11B = find(   p.f0mData(:,p.column_subject)==11 &    p.f0mData(:,p.column_SNR)==10);
idxS11C = find(p.aceDataOwn(:,p.column_subject)==11 & p.aceDataOwn(:,p.column_SNR)==10);
idxS13A = find(   p.aceData(:,p.column_subject)==13 &    p.aceData(:,p.column_SNR)==10);
idxS13B = find(   p.f0mData(:,p.column_subject)==13 &    p.f0mData(:,p.column_SNR)==10);
idxS13C = find(p.aceDataOwn(:,p.column_subject)==13 & p.aceDataOwn(:,p.column_SNR)==10);
idxS14A = find(   p.aceData(:,p.column_subject)==14 &    p.aceData(:,p.column_SNR)==10);
idxS14B = find(   p.f0mData(:,p.column_subject)==14 &    p.f0mData(:,p.column_SNR)==10);
idxS14C = find(p.aceDataOwn(:,p.column_subject)==14 & p.aceDataOwn(:,p.column_SNR)==10);

PlotMeans(stPlot, [p.aceData(idxS11A,p.column_score_phoneme)/MaxScore_p     p.aceData(idxS13A,p.column_score_phoneme)/MaxScore_p      p.aceData(idxS14A,p.column_score_phoneme   )/MaxScore_p]*100, ...
                  [p.f0mData(idxS11B,p.column_score_phoneme)/MaxScore_p     p.f0mData(idxS13B,p.column_score_phoneme)/MaxScore_p      p.f0mData(idxS14B,p.column_score_phoneme   )/MaxScore_p]*100);  

if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-noise-per-subject'],'epsc');
end

stPlot.bPlotIndividual = 0;
stPlot.xTick        = [1];
stPlot.xTickLabel   = {' '};

stPlot.TitleHead    = [' '];
stPlot.xLim         = [ 0   2];
stPlot.XLabel       = 'Pooled';

PlotMeans(stPlot,[p.aceData(idx_NoiseA ,p.column_score_phoneme   )/MaxScore_p]*100, ...
                 [p.f0mData(idx_NoiseB,p.column_score_phoneme   )/MaxScore_p]*100);  
stPlot.TitleHead    = ['']; 
stPlot.xTickLabel   = {' '};

if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-noise-pooled'],'epsc');
end



% Plotting LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjects = [11 13 14];
pLT = fullfil_subject_data(pLT, subjects);

idxS11A = find(pLT.aceData(:,pLT.column_subject)==11);
idxS11B = find(pLT.f0mData(:,pLT.column_subject)==11);
idxS11C = find(pLT.aceDataOwn(:,pLT.column_subject)==11);
idxS13A = find(pLT.aceData(:,pLT.column_subject)==13);
idxS13B = find(pLT.f0mData(:,pLT.column_subject)==13);
idxS13C = find(pLT.aceDataOwn(:,pLT.column_subject)==13);
idxS14A = find(pLT.aceData(:,pLT.column_subject)==14);
idxS14B = find(pLT.f0mData(:,pLT.column_subject)==14);
idxS14C = find(pLT.aceDataOwn(:,pLT.column_subject)==14);

stPlot.bPlotIndividual  = 1;
stPlot.xLim             = [ 0   3+2.5];
stPlot.yLim             = [-7   7];
stPlot.yTick            = [-7+1:2:7];
stPlot.xTick            = [1 2 3];
stPlot.xTickLabel       = {'RP','PM', 'JL'};
stPlot.XLabel           = 'Subject';
stPlot.YLabel           = 'SRT (dB)';

PlotMeans3series(stPlot,[pLT.aceData(idxS11A,2)    pLT.aceData(idxS13A,2)    pLT.aceData(idxS14A,2)], ...
                        [pLT.f0mData(idxS11B,2)    pLT.f0mData(idxS13B,2)    pLT.f0mData(idxS14B,2)], ...
                        [pLT.aceDataOwn(idxS11C,2) pLT.aceDataOwn(idxS13C,2) pLT.aceDataOwn(idxS14C,2)]);

if bSave
    saveas(gcf,[dest_folder,'LIST-per-subject'],'epsc');
end

stPlot.bPlotIndividual = 0;
stPlot.xTick        = [1];
stPlot.xTickLabel   = {' '};

stPlot.TitleHead    = [' '];
stPlot.xLim         = [ 0   2];
stPlot.XLabel       = 'Pooled';
PlotMeans3series(stPlot,[pLT.aceData(:,2)], ...
                        [pLT.f0mData(:,2)], ...
                        [pLT.aceDataOwn(:,2)]);  
stPlot.TitleHead    = ['']; 
stPlot.xTickLabel   = {' '};

if bSave
    saveas(gcf,[dest_folder,'LIST-pooled'],'epsc');
end

                    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MATRIX: 
% 
% stPlot.bPlotIndividual = 0;
% stPlot.TitleHead = ['Flemish Matrix, SNR = 10 dB (Results for CI 3 - Patrick Meul)']; 
% [handleFig pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);
% 
% if bSave
%     saveas(handleFig(1),[dest_folder,'MT-SNR10dB'],'epsc');
%     saveas(handleFig(2),[dest_folder,'MT-SNR05dB'],'epsc');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end