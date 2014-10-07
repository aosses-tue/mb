function experiment_report_20131023_roving
% function experiment_report_20131023_roving
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

hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
dest_folder     = [hoofd_folder 'ci_pooled/Figures/'];

SubjectName = 'CIs';

bPlot = 1;
bSave = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loudness Balancing

stPlot = [];
stPlot.SeriesLabel  = {'xPC F0mod', 'xPC ACE'};

% % When required type: 'ci-'
% [deltaACE, deltaF0mod]  = figLBScores_xPC(hoofd_folder,{'UW'}, bPlot);
% 
% if bSave
%     Min4plot = round(min(min([deltaACE; deltaF0mod]))-5);
%     Max4plot = round(max(max([deltaACE; deltaF0mod]))+5);
%     ylim([Min4plot Max4plot])
%     
%     saveas(gcf,[dest_folder,'LB-results'],'epsc');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pitch ranking

stSubject.firstname = 'ci-';
stSubject.truncate_data = 0;
stSubject.bRoving = 1;
stSubject.dest_folder = dest_folder;

[aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
% To confirm the data, type 1

if bSave
    saveas(gcf,[dest_folder,'PR-results-pooled'],'epsc');
end

stSubject.firstname = 'ci-Romain';
[aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
if bSave
    saveas(gcf,[dest_folder,'PR-results-RP'],'epsc');
end

stSubject.firstname = 'ci-Maria';
[aceData, f0mData]      = figPRScores_xPC(hoofd_folder,{'UW'}, stPlot, stSubject);
if bSave
    saveas(gcf,[dest_folder,'PR-results-MB'],'epsc');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end