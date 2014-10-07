function experiment_report_20130924
% function experiment_report_20130924
%
% Session with Romain Peeters, music experiments
%
% Dependencies:
%       figPRScores_xPC
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPlot = 1;
bSave = 0;

hoofd_folder        = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Romain_Peeters/';
dest_folder         = [hoofd_folder 'samenvatting-20130924/Figures/'];

result_folder       = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';

stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};

[deltaACE, deltaF0mod]  = figLBScores_xPC(result_folder,{'UW'}, bPlot);

if bSave
    close
    Min4plot = round(min(min([deltaACE; deltaF0mod]))-5);
    Max4plot = round(max(max([deltaACE; deltaF0mod]))+5);
    ylim([Min4plot Max4plot])
    xlim([100 300])
    saveas(gcf,[dest_folder,'LB-results'],'epsc');
end

[aceData, f0mData]      = figPRScores_xPC(result_folder,{'UW'}, stPlot);
% When required type: 'ci-Romain'
% To confirm the data, type 1

if bSave
    saveas(gcf,[dest_folder,'PR-results'],'epsc');
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end