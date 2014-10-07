function [stHandle] = quick_MT_one_subject(pMT, subject, nSessie)
% function [stHandle] = quick_MT_one_subject(pMT, subject, nSessie)
% 
% Get plots from Vlaamse Matrix (fixed SNRs at 10, 5 dB), for one of the 
% following subjects:
% 
% subject =     11 - Romain Peeters
%               12 - Maria Brughmans
%               13 - Patrick Meul
%               14 - Jan Leys
%               15 - Wouter David
%               16 - Jean Baptist Daumerie
%               17 - Julia Schoolmeesters
% % Example:
%       subject = 12;
%       stPlot.bPlotIndividual = 1;
%       file_MT_fixedSNR = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Maria_Brughmans/20131125-MT/20131125-CI-MB-VlMatrix.txt';
%       [xxx pMT] = figMTScores_xPC(file_MT_fixedSNR,stPlot);
%       quick_MT_one_subject(pMT, subject);
%
% Programmed by Alejandro Osses, ExpORL 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if subject < 10
    Subject = ['0' num2str(subject)];
else
    Subject = num2str(subject);
end

stPlot.TitleHead    = ['Flemish Matrix in noise, Subject S' Subject];
stPlot.xTickLabel   = {'+10','+5'};
stPlot.xTick        = 1:length(stPlot.xTickLabel);
stPlot.xLim         = [ 0   max(stPlot.xTick)];
stPlot.XLabel       = 'SNR (dB)';
stPlot.YLabel       = 'Word scores (%)';
stPlot.yLim         = [15 105];
stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
numDecimals         = 1; % Number of decimals to be displayed in figure legend
SNR = [10 5];
    
ace2plot        = [];
f0m2plot        = [];
ace2plot_tot    = [];
f0m2plot_tot    = [];

for k = 1:2

    tmpSSA          = [];
    tmpSSB          = [];

    idxSSA_tot  = find(   pMT.aceData(:,pMT.column_subject)==subject &    pMT.aceData(:,pMT.column_SNR)==SNR(k) );
    idxSSB_tot  = find(   pMT.f0mData(:,pMT.column_subject)==subject &    pMT.f0mData(:,pMT.column_SNR)==SNR(k) );
        
    if exist('nSessie','var')
        idxSSA      = find(   pMT.aceData(:,pMT.column_subject)==subject &    pMT.aceData(:,pMT.column_SNR)==SNR(k) & pMT.aceData(:,pMT.column_report)==nSessie);
        idxSSB      = find(   pMT.f0mData(:,pMT.column_subject)==subject &    pMT.f0mData(:,pMT.column_SNR)==SNR(k) & pMT.f0mData(:,pMT.column_report)==nSessie);
    else
        idxSSA = idxSSA_tot;
        idxSSB = idxSSB_tot;
    end
    
    ace2plot    = [ace2plot pMT.aceData(idxSSA,pMT.column_score_word)];
    f0m2plot    = [f0m2plot pMT.f0mData(idxSSB,pMT.column_score_word)];

    ace2plot_tot= [ace2plot_tot pMT.aceData(idxSSA_tot,pMT.column_score_word)];
    f0m2plot_tot= [f0m2plot_tot pMT.f0mData(idxSSB_tot,pMT.column_score_word)];
end
[mf0m sf0m]             = Get_mean(f0m2plot);
[mace sace]             = Get_mean(ace2plot);

mean2plot   = [mf0m' mace'];
std2plot    = [sf0m' sace'];

[mf0m_pooled sf0m_pooled] = Get_mean(f0m2plot_tot);
[mace_pooled sace_pooled] = Get_mean(ace2plot_tot);

stPlot.SeriesLabel  = { ['xPC F0mod'], ...
                         ['xPC ACE']};

mean2plot_pooled   = [mf0m_pooled' mace_pooled'];
std2plot_pooled    = [sf0m_pooled' sace_pooled'];

Colors = [0 0 0; 1 1 1];

figure
barweb(mean2plot,std2plot,[],[],stPlot.TitleHead,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel); hold on;

offset = 0.9/(2+1);
try 
    plot([1+offset/2 2+offset/2], ace2plot   ,'rx','LineWidth',2)
end
try
    plot([1-offset/2 2-offset/2], f0m2plot   ,'rx','LineWidth',2)
end

grid on
stHandle.h1 = gcf;
Handle = gca;
set(Handle,'YTick',stPlot.yTick)
set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
set(Handle,'XLim',[0 max(stPlot.xTick)+0.5])
set(Handle,'YLim',[10 110])
set(Handle,'XTickLabel',[stPlot.xTickLabel])

figure
barweb(mean2plot_pooled,std2plot_pooled,[],[],[stPlot.TitleHead ', pooled data'],stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
grid on
stHandle.hPooled = gcf;
Handle = gca;
set(Handle,'YTick',stPlot.yTick)
set(Handle,'XTick',[stPlot.xTick max(stPlot.xTick)+2])
set(Handle,'XLim',[0 max(stPlot.xTick)+0.5])
set(Handle,'YLim',[10 110])
set(Handle,'XTickLabel',[stPlot.xTickLabel])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end