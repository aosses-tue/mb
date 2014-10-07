function Matthias_experiments_dutch(dest_folder)
% function Matthias_experiments_dutch(dest_folder)

% Processes experiments conducted by Matthias during March-April 2009 (data unpublished).
%
% The experiments were:
%   - Sentence recognition in quiet
%   - Sentence recognition in noise (fixed SNR at 10 dB)
%
% Four CI subjects were tested, but only 3 were processed along this script
% (the information over the fourth subject was not found).
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDoSpeechTests  = 1;

if nargin == 0
    dest_folder = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan/proefschrift/Figures_new/';
end

result_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';

h = [];

if bDoSpeechTests
 
    [xx p] = figLTScores_xPC('CI-Sentence-LIST-fixed-MM.txt'); % 'CI-Sentence-recognition-MM.txt');
    
    if length(xx) ~= 0 % then figure was generated
        h(end+1) = xx;
    end
    % saveas(h(end),[dest_folder,'Sentence-results'],'epsc')

end

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_Quiet = find(p.aceData(:,2)==99);
idx_Noise = find(p.aceData(:,2)==10);

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
SeriesLabel         = {'NMT ACE', 'NMT F0mod'};
stPlot.TitleSuffix  = {''};
stPlot.SeriesLabel  = SeriesLabel;
stPlot.TitleHead    = ['LIST material (Results for 3 CI subject)']; 
stPlot.YLabel       = 'Sentence scores (%)';
stPlot.XLabel       = 'SNR (dB)';
stPlot.ReverseData  = 1; % from worst to best score
stPlot.separation_series = 0.1;

PlotMeans([p.aceData(idx_Quiet,1) p.aceData(idx_Noise,1)]/10*100, [p.f0mData(idx_Quiet,1) p.f0mData(idx_Noise,1)]/10*100, stPlot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end