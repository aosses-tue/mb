function experiment_report_20130926
% function experiment_report_20130926
%
% Session with Patrick Meul, speech tests
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

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Patrick_Meul/';
dest_folder     = [hoofd_folder 'samenvatting-20130926/Figures/'];

file_LIST_adaptive  = [hoofd_folder 'CI-Patrick-LT-adaptive.txt'];
file_VU_adaptive    = [hoofd_folder 'CI-Patrick-VU-adaptive.txt'];
file_MT_fixedSNR    = [hoofd_folder 'CI-Patrick-MT-fixedSNR.txt']
file_LL_fixedSNR    = [hoofd_folder 'CI-Patrick-LL-fixedSNR.txt'];

bPlot = 1;
bSave = 0;

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

stPlot.yLim = [0 10];
[xx pVUadaptive] = figVUScores_xPC(file_VU_adaptive, bPlot, stPlot);
if bSave
    saveas(gcf,[dest_folder,'VU-adaptive-results'],'epsc');
end

% Not possible to do analysis with 0 dof:
% [stat_p,stat_t,stat_st] = anova1([pVUadaptive.aceData(:,2) pVUadaptive.f0mData(:,2) pVUadaptive.aceDataOwn(:,2)]);
% [stat_c,stat_m,stat_h,stat_nms] = multcompare(stat_st);

[h p]       = figLLScores_xPC(file_LL_fixedSNR, stPlot);
close;

% Not possible to do analysis with 0 dof:
% % Without transformation:
% ToEv = [p.aceData(1:2,1) p.f0mData(1:2,1) p.aceDataOwn(:,1)]/10; % Words in quiet
% ToEv = [p.aceData(3:4,1) p.f0mData(3:4,1)]/10; % Words in noise
% ToEv = 2*asin(sqrt([p.aceData(1:2,1) p.f0mData(1:2,1) p.aceDataOwn(:,1)]/10)); % Word    in quiet
% ToEv = 2*asin(sqrt([p.aceData(1:2,2) p.f0mData(1:2,2) p.aceDataOwn(:,2)]/30)); % Phoneme in quiet
% ToEv = 2*asin(sqrt([p.aceData(1:2,3) p.f0mData(1:2,3) p.aceDataOwn(:,3)]/20)); % Consonant
% ToEv = 2*asin(sqrt([p.aceData(1:2,4) p.f0mData(1:2,4) p.aceDataOwn(:,4)]/10)); % Vowel
% ToEv = 2*asin(sqrt([p.aceData(3:4,1) p.f0mData(3:4,1)]/10)); % Word    in quiet
% ToEv = 2*asin(sqrt([p.aceData(3:4,2) p.f0mData(3:4,2)]/30)); % Phoneme in quiet
% ToEv = 2*asin(sqrt([p.aceData(3:4,3) p.f0mData(3:4,3)]/20)); % Consonant
% ToEv = 2*asin(sqrt([p.aceData(3:4,4) p.f0mData(3:4,4)]/10)); % Vowel
% [stat_p,stat_t,stat_st] = anova1(ToEv);
% [stat_c,stat_m,stat_h,stat_nms] = multcompare(stat_st);


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
stPlot.yLim         = [35 105];
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

stPlot.TitleHead    = ['Lilliput in quiet (Results for CI 3 - Patrick Meul)']; 
stPlot.figPos       = [0 0 1024 450];
stPlot.xLim         = [ 0   5];
stPlot.xTick        = [1 2 3 4];
stPlot.xTickLabel   = {'Word','Phoneme','Vowel', 'Cons.'};
stPlot.XLabel       = 'Score type';
stPlot.YLabel       = 'Scores (%)';
score2plot          = p.column_score_vowel;
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

stPlot.TitleHead    = ['Lilliput, SNR = 10 dB (Results for CI 3 - Patrick Meul)']; 
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
% MATRIX: 

stPlot.bPlotIndividual = 0;
stPlot.TitleHead = ['Flemish Matrix, SNR = 10 dB (Results for CI 3 - Patrick Meul)']; 
[handleFig pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);

if bSave
    saveas(handleFig(1),[dest_folder,'MT-SNR10dB'],'epsc');
    saveas(handleFig(2),[dest_folder,'MT-SNR05dB'],'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end