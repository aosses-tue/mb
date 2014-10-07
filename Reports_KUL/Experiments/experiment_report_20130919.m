function experiment_report_20130919
% function experiment_report_20130919
%
% Session with Romain Peeters
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

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Romain_Peeters/';
dest_folder     = [hoofd_folder 'samenvatting-20130919/Figures/'];

file_LIST_adaptive  = [hoofd_folder 'CI-Romain-LT-adaptive.txt'];
file_VU_adaptive    = [hoofd_folder 'CI-Romain-VU-adaptive.txt'];
file_VU_fixedSNR    = [hoofd_folder 'CI-Romain-VU-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'CI-Romain-Word-LL-fixedSNR.txt'];

bPlot = 1;
bSave = 1;

stPlot.xTick            = 1;
stPlot.bPlotIndividual  = 1;
stPlot.YGrid            = 'on';
stPlot.yLim             = [-6 4];
[xx pLT]                = figLTScores_xPC(file_LIST_adaptive, bPlot, stPlot);

if bSave
    saveas(gcf,[dest_folder,'LT-adaptive-results'],'epsc');
end

% [stat_p,stat_t,stat_st] = anova1([pLT.aceData(:,2) pLT.f0mData(:,2) pLT.aceDataOwn(:,2)]);
% [stat_c,stat_m,stat_h,stat_nms] = multcompare(stat_st);
% Result: no groups have signifficantly different (p=73,16%)

stPlot.yLim = [0 20];
[xx pVUadaptive] = figVUScores_xPC(file_VU_adaptive, bPlot, stPlot);
if bSave
    saveas(gcf,[dest_folder,'VU-adaptive-results'],'epsc');
end

% [stat_p,stat_t,stat_st] = anova1([pVUadaptive.aceData(:,2) pVUadaptive.f0mData(:,2) pVUadaptive.aceDataOwn(:,2)]);
% [stat_c,stat_m,stat_h,stat_nms] = multcompare(stat_st);
% Result: no groups have signifficantly different (p=%)

[xx pVU]    = figSentenceScores_xPC(file_VU_fixedSNR);

[h p]       = figLLScores_xPC(file_LL_fixedSNR, stPlot);
close;

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_Quiet = find(p.aceData(:,p.column_SNR)==99);
idx_Noise = find(p.aceData(:,p.column_SNR)==10);

stPlot.nCols        = 1;
stPlot.nRows        = 1;
stPlot.xGap         = 15;
stPlot.yGap         = 80;
stPlot.fntsz        = 16;
stPlot.markerSize   = 12;
stPlot.margins      = [70 70 50 90];
stPlot.figPos       = [0 0 1024 450]; % 225 each row % left, right, top, bottom
stPlot.xLim         = [  0  3];
stPlot.yLim         = [0 100];
stPlot.yTick        = [stPlot.yLim(1):10:stPlot.yLim(2)];
stPlot.xTickLabel   = {'Q','+10'};
stPlot.YGrid        = 'on';
stPlot.TitleSuffix  = {''};
stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod', 'ACE baseline'};
stPlot.XLabel       = 'SNR (dB)';
stPlot.ReverseData  = 0; % from worst to best score
stPlot.separation_series = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lilliput scores: 

stPlot.TitleHead    = ['Lilliput in quiet (Results for CI 1 - Romain Peeters)']; 
stPlot.figPos       = [0 0 1024 450];
stPlot.xLim         = [ 0   5];
stPlot.xTick        = [1 2 3 4];
stPlot.xTickLabel   = {'Word','Phoneme','Vowel', 'Cons.'};
stPlot.XLabel       = 'Score type';
stPlot.YLabel       = 'Scores (%)';
%score2plot          = p.column_score_vowel;
MaxScore_v          = 10;
MaxScore_c          = 20;
MaxScore_p          = 30;
MaxScore_w          = 10;

PlotMeans3series(stPlot,[p.aceData(idx_Quiet,p.column_score_word   )/MaxScore_w   p.aceData(idx_Quiet,p.column_score_phoneme   )/MaxScore_p   p.aceData(idx_Quiet,p.column_score_vowel   )/MaxScore_v  p.aceData(idx_Quiet,p.column_score_cons)/MaxScore_c]*100, ...
                        [p.f0mData(idx_Quiet,p.column_score_word   )/MaxScore_w   p.f0mData(idx_Quiet,p.column_score_phoneme   )/MaxScore_p   p.f0mData(idx_Quiet,p.column_score_vowel   )/MaxScore_v  p.f0mData(idx_Quiet,p.column_score_cons)/MaxScore_c]*100, ...
                        [p.aceDataOwn(idx_Quiet,p.column_score_word)/MaxScore_w   p.aceDataOwn(idx_Quiet,p.column_score_phoneme)/MaxScore_p   p.aceDataOwn(idx_Quiet,p.column_score_vowel)/MaxScore_v  p.aceDataOwn(idx_Quiet,p.column_score_cons)/MaxScore_c]*100);  
if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-quiet'],'epsc');
end

stPlot.TitleHead    = ['Lilliput, SNR = 10 dB (Results for CI 1 - Romain Peeters)']; 
stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};
stPlot.xTickLabel   = {'Word', 'Phoneme','Vowel', 'Cons.'};
stPlot.xLim         = [  0  5];
stPlot.line_color   = 'w'; % no line
stPlot.bPlotIndividual = 1;

PlotMeans(stPlot, [p.aceData(idx_Noise,p.column_score_word)/MaxScore_w  p.aceData(idx_Noise,p.column_score_phoneme)/MaxScore_p   p.aceData(idx_Noise,p.column_score_vowel   )/MaxScore_v  p.aceData(idx_Noise,p.column_score_cons)/MaxScore_c]*100, ...
                  [p.f0mData(idx_Noise,p.column_score_word)/MaxScore_w  p.f0mData(idx_Noise,p.column_score_phoneme)/MaxScore_p   p.f0mData(idx_Noise,p.column_score_vowel   )/MaxScore_v  p.f0mData(idx_Noise,p.column_score_cons)/MaxScore_c]*100);  

if bSave
    saveas(gcf,[dest_folder,'Lilliput-p-v-c-scores-noise'],'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VU Sentence scores: 

idx_SNR05           = find(pVU.aceData(:,pVU.column_SNR)== 5);
idx_SNR10           = find(pVU.aceData(:,pVU.column_SNR)==10);

stPlot.TitleHead    = ['VU sentences, SNR = 10 dB (Results for CI 1 - Romain Peeters)']; 
stPlot.xLim         = [ 0   3];
stPlot.xTick        = [1 2];
stPlot.xTickLabel   = {'Sentence','Word'};
stPlot.XLabel       = 'Score type';
stPlot.YLabel       = 'Scores (%)';
MaxScore_s          = 13;
MaxScore_w          = 1; % different lists with different number of words
stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};

PlotMeans(stPlot,[pVU.aceData(idx_SNR10,pVU.column_score_sentence)/MaxScore_s   pVU.aceData(idx_SNR10,pVU.column_score_word)/MaxScore_w]*100, ...
                 [pVU.f0mData(idx_SNR10,pVU.column_score_sentence)/MaxScore_s   pVU.f0mData(idx_SNR10,pVU.column_score_word)/MaxScore_w]*100)
close
if bSave
    %saveas(gcf,[dest_folder,'VU-SNR10'],'epsc');
end

stPlot.TitleHead    = ['VU sentences in noise (Results for CI 1 - Romain Peeters)']; 
PlotMeans(stPlot,[pVU.aceData(idx_SNR05,pVU.column_score_sentence)/MaxScore_s   pVU.aceData(idx_SNR05,pVU.column_score_word)/MaxScore_w]*100, ...
                 [pVU.f0mData(idx_SNR05,pVU.column_score_sentence)/MaxScore_s   pVU.f0mData(idx_SNR05,pVU.column_score_word)/MaxScore_w]*100)
if bSave
    %saveas(gcf,[dest_folder,'VU-SNR05'],'epsc');
end
close
stPlot.xLim         = [ 0   5.5];
stPlot.xTick        = [1 2 3 4];
stPlot.xTickLabel   = {'Sentence, +10','Sentence, +5','Word, +10','Word, +5'};
stPlot.XLabel       = 'Score type, SNR (dB)';
stPlot.bPlotIndividual = 0;

PlotMeans(stPlot,[pVU.aceData(idx_SNR10,pVU.column_score_sentence)/MaxScore_s pVU.aceData(idx_SNR05,pVU.column_score_sentence)/MaxScore_s   pVU.aceData(idx_SNR10,pVU.column_score_word)/MaxScore_w pVU.aceData(idx_SNR05,pVU.column_score_word)/MaxScore_w]*100, ...
                 [pVU.f0mData(idx_SNR10,pVU.column_score_sentence)/MaxScore_s pVU.f0mData(idx_SNR05,pVU.column_score_sentence)/MaxScore_s   pVU.f0mData(idx_SNR10,pVU.column_score_word)/MaxScore_w pVU.f0mData(idx_SNR05,pVU.column_score_word)/MaxScore_w]*100)
hold on
plot([2.5 2.5], [0 100],'k')

if bSave
    saveas(gcf,[dest_folder,'VU-fixed-SNR'],'epsc');
end

% Plotting Staircases VU:

VU_List_18 = [0 2 4 6 8 10 12 14 12 14 12 14 12 10 12 10 12 10 12 10 12];
quick_Sentence_adaptive(VU_List_18,1);
ylim([-2 22])
handle = gca;
set(handle,'YTick',-2:2:22)

if bSave
    saveas(gcf,[dest_folder,'Staircase-F0mod-1'],'epsc');
end

VU_List_19 = [0 2 4 6 8 10 8 6 8 6 8 10 8 10 12 14 16 18 20];
quick_Sentence_adaptive(VU_List_19,1);
handle = gca;
set(handle,'YTick',-2:2:22)

if bSave
    saveas(gcf,[dest_folder,'Staircase-F0mod-2'],'epsc');
end

VU_List_22 = [0 2 4 6 8 10 12 10 8 10 12 10 8 6 8 6 4 2 4 6];
quick_Sentence_adaptive(VU_List_22,0);

if bSave
    saveas(gcf,[dest_folder,'Staircase-ACE-1'],'epsc');
end

VU_List_23 = [0 2 4 6 8 10 12 10 8 6 8 6 4 6 8 10 8 6 8 10];
quick_Sentence_adaptive(VU_List_23,0);
 
if bSave
    saveas(gcf,[dest_folder,'Staircase-ACE-2'],'epsc');
end

VU_List_16 = [0 2 4 6 4 2 4 6 4 6 8 6 8 10 12 10 8];
quick_Sentence_adaptive(VU_List_16,10);
 
if bSave
    saveas(gcf,[dest_folder,'Staircase-ACEOwn-1'],'epsc');
end

VU_List_17 = [0 2 4 2 4 6 8 10 12 10 8 10 12 10 8 6];
quick_Sentence_adaptive(VU_List_17,10);
 
if bSave
    saveas(gcf,[dest_folder,'Staircase-ACEOwn-2'],'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end