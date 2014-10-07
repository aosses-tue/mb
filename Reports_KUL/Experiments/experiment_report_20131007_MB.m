function experiment_report_20131007_MB
% function experiment_report_20131007_MB
%
% Session with Maria Brughmans, music experiments
%
% Dependencies:
%       figPRScores_xPC
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPlot = 1;
bSave = 1;

hoofd_folder        = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Maria_Brughmans/';
dest_folder         = [hoofd_folder 'samenvatting-20131007/Figures/'];

result_folder       = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';

stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};

[deltaACE, deltaF0mod]  = figLBScores_xPC(result_folder,{'UW'}, bPlot);
close

if bSave
    close
    Min4plot = round(min(min([deltaACE; deltaF0mod]))-5);
    Max4plot = round(max(max([deltaACE; deltaF0mod]))+5);
    ylim([Min4plot Max4plot])
    xlim([100 300])
    saveas(gcf,[dest_folder,'LB-results'],'epsc');
end

[aceData, f0mData]      = figPRScores_xPC(result_folder,{'UW'}, stPlot);
% When required type: 'ci-Maria'
% To confirm the data, type 1

if bSave
    saveas(gcf,[dest_folder,'PR-results'],'epsc');
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end