function [h stPlot] = PlotMeans_default(stPlot, aceData, f0mData, aceData2, f0mData2)
% function [h stPlot] = PlotMeans_default(stPlot, aceData, f0mData)
%
%
% stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};
% stPlot.TitleHead    = ['LIST material']; 
% stPlot.YLabel       = 'SRT (dB)';
% stPlot.xLim     = [0  size(aceData,2)];
% stPlot.xLim     = [0 size(aceData,2)*4+1];
% stPlot.xTick    = min(stPlot.xLim)+1:1:max(stPlot.xLim)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    stPlot = [];
end
h = [];

stPlot = Ensure_field(stPlot,'nCols',1);
stPlot = Ensure_field(stPlot,'nRows',1);
stPlot = Ensure_field(stPlot,'xGap' ,15);
stPlot.yGap         = 80;
stPlot.fntsz        = 16;
stPlot.markerSize   = 12;
stPlot.margins      = [70 70 50 90];
stPlot.figPos       = [0 0 1024 450]; % 225 each row % left, right, top, bottom
if nargin >= 3
    stPlot.xLim     = [0  size(aceData,2)];
    stPlot.xLim     = [0 size(aceData,2)*4+1];
    stPlot.xTick    = min(stPlot.xLim)+1:1:max(stPlot.xLim)-1;
end

stPlot.yLim         = [-10 10];
stPlot.yTick        = -8:2:8;
stPlot.xTickLabel   = {''};
stPlot.YGrid        = 'on';
stPlot.TitleSuffix  = {''};
stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod'};
stPlot.TitleHead    = ['LIST material']; 
stPlot.YLabel       = 'SRT (dB)';
stPlot.XLabel       = ' ';
% stPlot.ReverseData  = 1; % from worst to best score
stPlot.separation_series = 0.1;

if nargin == 3
    PlotMeans(aceData, f0mData, stPlot);
elseif nargin > 3
    try
        PlotMeans(aceData, f0mData, stPlot, aceData2, f0mData2);
    catch
        stPlot.xTickLabel   = {'1. ACE','1. F0mod, IIR ED','2. ACE', '2. F0mod, FIR ED'};
        stPlot.SeriesLabel  = {'xPC ACE', 'xPC F0mod','minmax SRT'};
        Avg        = [mean([aceData f0mData]) mean([aceData2 f0mData2])];
        Std        = [std([aceData f0mData]) std([aceData2 f0mData2])];
        AvgACE     = [mean(aceData) mean(aceData2)];
        StdACE     = [std(aceData)  std(aceData2)];
        AvgF0m     = [mean(f0mData) mean(f0mData2)];
        StdF0m     = [std(f0mData)  std(f0mData2)];
        h = figure;
        set(h,'Position',stPlot.figPos);

        hp = plot([1 3], AvgACE,'Marker', 'o','Color',[0.75 0.75 0.75],'LineWidth',4,'LineStyle','none'); hold on, grid on
        plot([2 4], AvgF0m, 'ko','LineWidth',4)

        ha      = gca;
        set(ha,'XTick'      , stPlot.xTick      );
        set(ha,'XTickLabel' , stPlot.xTickLabel );
        set(ha,'YTick'      , stPlot.yTick      );

        xlabel( stPlot.XLabel );
        ylabel( stPlot.YLabel );

        title ( ['Speech-in-noise scores using ' stPlot.TitleHead ' (4 Subjects)'] );
        plot([minmax(aceData'); minmax(f0mData'); minmax(aceData2'); minmax(f0mData2')] ,'rx','LineWidth',2)
        legend(stPlot.SeriesLabel)

        errorbar(Avg,Std,'k.'), hold on
        xlim(stPlot.xLim)
        ylim(stPlot.yLim)
    end
elseif nargin <= 1
    disp([mfilename '.m: you may need to define more parameters to plot. When ready use the command ''PlotMeans(data1,  data2, stPlot)'''])
end

if nargin >= 3
    set(gcf, 'PaperPositionMode','auto')
end

end