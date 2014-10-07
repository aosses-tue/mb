function [stHandle] = quick_LL_one_subject(p, subject, nSessie)
% function [stHandle] = quick_LL_one_subject(p, subject, nSessie)
% 
% Get plots from Lilliput (quiet and SNR at 10 dB), for one of the 
% following subjects:
% 
% subject =     11 - Romain Peeters
%               12 - Maria Brughmans
%               13 - Patrick Meul
%               14 - Jan Leys
%               15 - Wouter David
%               16 - Jean Baptist Daumerie
%               17 - Julia Schoolmeesters
% 
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if subject < 10
    Subject = ['0' num2str(subject)];
else
    Subject = num2str(subject);
end

% score2use           = p.column_score_word;
score2use           = p.column_score_phoneme;
value2norm          = 30/100; % out of 30 and express in percentage

stPlot.TitleHead    = ['Lilliput, Subject S' Subject];
stPlot.xTickLabel   = {'Q','+10'};
stPlot.xTick        = 1:length(stPlot.xTickLabel);
stPlot.xLim         = [ 0   max(stPlot.xTick)];
stPlot.XLabel       = 'SNR (dB)';
stPlot.YLabel       = 'Phoneme scores (%)';
stPlot.yLim         = [15 105];
stPlot.yTick        = [stPlot.yLim(1)+5:10:stPlot.yLim(2)];
SNR                 = [99 10];

ace2plot        = [];
f0m2plot        = [];
aceOwn2plot     = [];
ace2plot_tot    = [];
f0m2plot_tot    = [];
aceOwn2plot_tot = [];

for k = 1:2

    idxSSA_tot  = find(   p.aceData(:,p.column_subject)==subject &    p.aceData(:,p.column_SNR)==SNR(k) );
    idxSSB_tot  = find(   p.f0mData(:,p.column_subject)==subject &    p.f0mData(:,p.column_SNR)==SNR(k) );
    idxSSC_tot  = find(   p.f0mData(:,p.column_subject)==subject &    p.f0mData(:,p.column_SNR)==SNR(k) );
        
    if exist('nSessie','var')
        idxSSA      = find(   p.aceData(:,p.column_subject)==subject &    p.aceData(:,p.column_SNR)==SNR(k) & p.aceData(:,p.column_report)==nSessie);
        idxSSB      = find(   p.f0mData(:,p.column_subject)==subject &    p.f0mData(:,p.column_SNR)==SNR(k) & p.f0mData(:,p.column_report)==nSessie);
        idxSSC      = find(   p.aceDataOwn(:,p.column_subject)==subject &    p.aceDataOwn(:,p.column_SNR)==SNR(k) & p.aceDataOwn(:,p.column_report)==nSessie);
    else
        idxSSA = idxSSA_tot;
        idxSSB = idxSSB_tot;
        idxSSC = idxSSC_tot;
    end
    
    try
        ace2plot        = [ace2plot p.aceData(idxSSA,score2use)/value2norm];
    catch
        if length(idxSSA) > length(ace2plot)
            ace2plot    = [ace2plot; nan(length(idxSSA)-length(ace2plot),1)];
            ace2plot    = [ace2plot p.aceData(idxSSA,score2use)/value2norm];
        end
    end
    
    try
        f0m2plot        = [f0m2plot p.f0mData(idxSSB,score2use)/value2norm];
    catch
        if length(idxSSB) > length(f0m2plot)
            f0m2plot    = [f0m2plot; nan(length(idxSSB)-length(f0m2plot),1)];
            f0m2plot    = [f0m2plot p.f0mData(idxSSB,score2use)/value2norm];
        end
    end
    
    try
        aceOwn2plot    = [aceOwn2plot p.aceDataOwn(idxSSC,score2use)/value2norm];
    catch 
        if length(idxSSC) > length(aceOwn2plot)
            aceOwn2plot = [aceOwn2plot; nan(length(idxSSC)-length(aceOwn2plot),1)];
            aceOwn2plot    = [aceOwn2plot p.aceDataOwn(idxSSC,score2use)/value2norm];
        end
    end

    ace2plot_tot= [ace2plot_tot p.aceData(idxSSA_tot,score2use)/value2norm];
    f0m2plot_tot= [f0m2plot_tot p.f0mData(idxSSB_tot,score2use)/value2norm];
    aceOwn2plot_tot= [aceOwn2plot_tot p.aceDataOwn(idxSSC_tot,score2use)/value2norm];
end
[mf0m sf0m]             = Get_mean(f0m2plot);
[mace sace]             = Get_mean(ace2plot);
[maceOwn saceOwn]       = Get_mean(aceOwn2plot);

mean2plot   = [mf0m' mace' maceOwn'];
std2plot    = [sf0m' sace' saceOwn'];

[mf0m_pooled sf0m_pooled] = Get_mean(f0m2plot_tot);
[mace_pooled sace_pooled] = Get_mean(ace2plot_tot);
[maceOwn_pooled saceOwn_pooled] = Get_mean(aceOwn2plot_tot);

stPlot.SeriesLabel  = { 'xPC F0mod', 'xPC ACE', 'baseline'};

mean2plot_pooled   = [mf0m_pooled' mace_pooled' maceOwn_pooled'];
std2plot_pooled    = [sf0m_pooled' sace_pooled' saceOwn_pooled'];

Colors = [0 0 0; 1 1 1; 0.75 0.75 0.75];

figure
barweb(mean2plot,std2plot,[],[],stPlot.TitleHead,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel), hold on;

offset = 0.9/(3+1);
try
    plot([1        2       ], ace2plot   ,'rx','LineWidth',2)
end
try
    plot([1-offset 2-offset], f0m2plot   ,'rx','LineWidth',2)
end
try
    plot([1+offset 2+offset], aceOwn2plot,'rx','LineWidth',2)
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