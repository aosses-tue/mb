function ac = PlotMeans3series(stPlot, x1, x2, x3) 
% function PlotMeans3series(stPLot, x1, x2, x3)
%
% x1, x2 and the optional x3 and x4 should have the same size. Each row will 
% correspond to a different subject.
%
% stPlot is a struct containing the Figure format:
%   xTicks - numeric values you want to display along the x-axis, its length
%            defines the number of points to be plotted per serie
%
% Recommended fields to be defined:
%
% nRepetitions        = p.nRepetitions; 
% nStrategies         = p.nStrategies;
% nSemitones          = p.nSemitones;
%
% stPlot.nPlots       = nReferences;
% stPlot.nCols        = 2;
% stPlot.nRows        = 1;
% stPlot.xGap         = 15;
% stPlot.yGap         = 80;
% stPlot.fntsz        = 16;
% stPlot.markerSize   = 12;
% stPlot.margins      = [70 70 50 90];
% stPlot.figPos       = [0 0 1024 225*nRows];
% stPlot.xTick        = 1:nSemitones;
% stPlot.yTick        = 0:20:100;
% [tones, xTicks]     = get_tones_PR(regs);   stPlot.xTickLabel       = xTicks(:,2:end);
% slev    = get_significance_level(nSubjects*nRepetitions, nStrategies); stPlot.slev = slev;
% stPlot.SeriesLabel  = {'ACE', 'F0mod'};
% TitleHead           = ['Pitch Ranking, Ref. '];     stPlot.TitleHead = TitleHead;
% stPlot.YLabel       = '% Correct';
% stPlot.XLabel      1 = 'Comparison tones [Hz]';
% stPlot.ReverseData  = 1; % from worst to best score
% stPlot.xTick        = xTick;
% stPlot.yTick        = yTick;
% stPlot.TitleSuffix  = {'PR ref 131 [Hz]',''PR ref 208 [Hz]''};
%
% stPlot.bTruncateName = 0;
%
% PlotMeans(aceData, f0mData, stPlot);
%
% Examples: 
%       figPRScores_xPC('/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/',{'CLARINET'});
%
% Plots using the Matthias' MATLAB classes
%
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    stPlot = [];
end

stPlot = Ensure_field(stPlot,'xTick',1:size(x1,2));
stPlot = Ensure_field(stPlot,'yTick',0:20:100);
stPlot = Ensure_field(stPlot,'nPlots', 1);
stPlot = Ensure_field(stPlot,'scalerate',1); % No scaling

try
    nReferences = length(stPlot.TitleSuffix);
catch
    disp('Using default number of Reference tones')
    nReferences  = 4; % Gives the number of plots
end

if nargin <= 3
    stPlot = Ensure_field(stPlot,'ReverseData', 1);
end

stPlot = Ensure_field(stPlot,'nCols'         , 1);
stPlot = Ensure_field(stPlot,'nRows'         , ceil(stPlot.nPlots/stPlot.nCols)); 
stPlot = Ensure_field(stPlot,'xGap'          , 10);
stPlot = Ensure_field(stPlot,'yGap'          , 60);
stPlot = Ensure_field(stPlot,'margins'       , [70 70 50 70]);
stPlot = Ensure_field(stPlot,'figPos'        , [0 0 1024 768/2]);

switch nargin
    case 4
        stPlot = Ensure_field(stPlot,'SeriesColor'   , {'k'  , 'w', [0.8 0.8 0.8]});
    otherwise
        stPlot = Ensure_field(stPlot,'SeriesColor'   , {'k'  , 'w'});
end
    
stPlot = Ensure_field(stPlot,'SeriesLabel'   , {'Serie 1','Serie 2','Serie 3'});
stPlot = Ensure_field(stPlot,'TitleHead'     , ['']);
stPlot = Ensure_field(stPlot,'TitleSuffix'   , ['']);
stPlot = Ensure_field(stPlot,'XLabel'        , ['Default x-axis name']);
stPlot = Ensure_field(stPlot,'YLabel'        , ['Default y-axis name']);
stPlot = Ensure_field(stPlot,'xTickLabel'    , stPlot.xTick);
stPlot = Ensure_field(stPlot,'yTickLabel'    , stPlot.yTick);
stPlot = Ensure_field(stPlot,'sameYAxis'     , 1);
stPlot = Ensure_field(stPlot,'LocationLegend', 'SouthEast');
stPlot = Ensure_field(stPlot,'ReverseData'   , 0);
stPlot = Ensure_field(stPlot,'separation_series', 0.1);
stPlot = Ensure_field(stPlot,'XGrid'         ,'off');
stPlot = Ensure_field(stPlot,'YGrid'         ,'on');
stPlot = Ensure_field(stPlot,'bPlotIndividual',0);
if stPlot.sameYAxis == 0
    stPlot.xGap = 110;
end
stPlot = Ensure_field(stPlot,'bTruncateName' , 0);

ColorLines  = 'k'; 
fntsz       = 12;

CountFig    = 1:stPlot.nPlots;
nRows       = stPlot.nRows;
nCols       = stPlot.nCols;
xGap        = stPlot.xGap;
yGap        = stPlot.yGap;
margins     = stPlot.margins;
figPos      = stPlot.figPos;
xTickLabel  = stPlot.xTickLabel;
yTickLabel  = stPlot.yTickLabel;
xTick       = stPlot.xTick;
yTick       = stPlot.yTick;

separation_series = stPlot.separation_series;

if length(stPlot.TitleSuffix) == 1
    try
        TitleSuffix = name2figname( stPlot.TitleSuffix );
    catch
        TitleSuffix = stPlot.TitleSuffix;
    end
else
    for i = 1:length(stPlot.TitleSuffix)
        TitleSuffix{i} = name2figname( stPlot.TitleSuffix{i} );
    end
end

SeriesColor = stPlot.SeriesColor;
SeriesLabel = stPlot.SeriesLabel;
XLabel      = stPlot.XLabel;
YLabel      = stPlot.YLabel;
TitleHead   = name2figname( stPlot.TitleHead );
LocationLegend = stPlot.LocationLegend;
markerSize   = round(0.005*figPos(3)); % A multiple of the window height %12;

pc  = setupPlotConf(nRows, nCols, false, xGap, yGap, margins, figPos, stPlot.YGrid);

stPlot = Ensure_field(stPlot,'xLim',[min(xTick)-0.2 max(xTick)+0.2]);
try
    stPlot = Ensure_field(stPlot,'yLim',[min(yTick)-20  max(yTick)+5]);
catch
    stPlot = Ensure_field(stPlot,'yLim',[min(yTick{1,1})-20  max(yTick{1,1})+5]);
    display('yLim: maybe stPlot.yTick is not a double variable, but a cell')
end

xLim    = stPlot.xLim;
yLim    = stPlot.yLim;

if nCols >= 2
    if mod(nReferences, 2) == 1
        set(pc.axesConfig(nRows,2).handle, 'visible', 'off')
    end
end
    

CountRow = [1 1 2 2 3 3]; % Horizontal layer
CountCol = [1 2 1 2 1 1]; % Vertical layer
CountID  = {'A','B','C','D','E','F'}; % No more than 6 plots per page

for i = CountFig
    
    if stPlot.ReverseData == 1
        try
            DataSerie1  = x1(:, (length(xTick):-1:1) + length(xTick)*(i-1) ); % Reverse the data (before: from interval 4 to 1, after: from 1 to 4)
            DataSerie2  = x2(:, (length(xTick):-1:1) + length(xTick)*(i-1) );
            if nargin == 4
                DataSerie3 = x3(:, (length(xTick):-1:1) + length(xTick)*(i-1) ); 
            end
        catch
            error([mfilename ': mismatch in data length, check that x1 and x2 have the same length than stPlot.xTick'])
        end
    else
        DataSerie1 = x1;
        DataSerie2 = x2;
        DataSerie3 = x3;
    end
    
    %%% Up to here depends on the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = CountRow(i);
    col = CountCol(i);
    ID  = CountID{i};
    
    ac  = pc.axesConfig(row,col);
    
    try
        ac.ax_xLim = xLim(i,:);
    catch
        ac.ax_xLim = xLim;
    end
    
    try
        ac.ax_yLim = yLim(i,:);  
    catch
        ac.ax_yLim = yLim;
    end 
    
    if size(xTickLabel,1) == 1 % only one row
        ac.ax_xTickLabel = xTickLabel(1,:);  % This could be automated
    else
        ac.ax_xTickLabel = xTickLabel(i,:);  % This could be automated
    end
    
    try
        ac.ax_xTick = xTick;
    catch
        ac.ax_xTick = xTick{1};
    end

    try
        ac.ax_yTick = yTick;
    catch
        ac.ax_yTick = yTick{1};
    end

    if mod(i,2) ~= 0 % then odd number
        if size(yTickLabel,1) == 1 % only one row
            ac.ax_yTickLabel = yTickLabel(1,:);  % This could be automated
        else
            ac.ax_yTickLabel = yTickLabel(i,:);  % This could be automated
        end
    end
    
    pcplt1 = PercentCorrectPlotter(     ac, ...                            % ac is the figure object
                                        DataSerie1, ... % aceData (size = N x M), N subjects, M different trials
                                        -2*separation_series, ...
                                        0, ...
                                        {}); 
    pcplt1.setupMarkers('s', markerSize, SeriesColor{1}, SeriesLabel{1}); % 1 = ACE
    pcplt1.setupLines(':', 1, 'w', '');
    pcplt1.setupErrorbars('-', 1, ColorLines, 0.1);
    pcplt1.plot();
    
    if nargin > 2
    pcplt2 = PercentCorrectPlotter(     ac, ...
                                        DataSerie2, ...
                                        0, ...
                                        0, ...
                                        {});
    pcplt2.setupMarkers('s', markerSize, SeriesColor{2}, SeriesLabel{2});
    pcplt2.setupLines(':', 1, 'w', '');
    pcplt2.setupErrorbars('-', 1, ColorLines, 0.1);
    pcplt2.plot()
    end
    
    if nargin > 3
    pcplt3 = PercentCorrectPlotter(     ac, ...
                                        DataSerie3, ...
                                        2*separation_series, ...
                                        0, ...
                                        {});
    pcplt3.setupMarkers('s', markerSize, SeriesColor{3}, SeriesLabel{3});
    pcplt3.setupLines(':', 1, 'w', '');
    pcplt3.setupErrorbars('-', 1, ColorLines, 0.1);
    pcplt3.plot()
    end
    
    if stPlot.bPlotIndividual == 1
        hold on;
        xTemp = ( 1:size(DataSerie1,2) )-2*separation_series;
        plot(xTemp, DataSerie1, 'rx', 'LineWidth',0.2)
        xTemp = ( 1:size(DataSerie2,2) )-0;
        plot(xTemp, DataSerie2, 'rx', 'LineWidth',0.2)
        xTemp = ( 1:size(DataSerie3,2) )+2*separation_series;
        try
            plot(xTemp, DataSerie3, 'rx', 'LineWidth',0.2)
        end
        if strcmp(stPlot.YGrid,'on')
            set(gca,'YGrid','on');
        end
    end
    
    ac.ax_box = 'on';
    
    h1 = findobj(ac.handle, 'tag', SeriesLabel{1}); % 1 = ACE
    h2 = findobj(ac.handle, 'tag', SeriesLabel{2}); % 2 = F0m
    h3 = findobj(ac.handle, 'tag', SeriesLabel{3});
    
    label = []; % Refresh label for axis
    positionLabel = 'middle';
    
    if stPlot.sameYAxis == 0 || i==1
        legend([h1(1) h2(1) h3(1)], SeriesLabel, 'location', LocationLegend);
        try
            yl1 = YLabelConfig(pc.axesConfig(1,i), YLabel{i}, yGap+60, 'left', 'middle', 90); 
        catch % Then YLabel is not a cell array
            yl1 = YLabelConfig(pc.axesConfig(1,1), YLabel, yGap, 'left', 'middle', 90); 
        end
    elseif mod(i,2) == 1
        yl1 = YLabelConfig(pc.axesConfig(row,1), YLabel, yGap, 'left', 'middle', 90); 
    end
    yl1.plot();
    
    if i == 1
        label = TitleHead;
    end
    
    try
        label = [label ' ' stPlot.TitleSuffix{i}];
    catch
        label = [label ' ' stPlot.TitleSuffix];
    end
    
    if stPlot.bTruncateName == 1
        if length(label)>40
            label = label(1:40);
            display('Filename truncated while plotting')
        end
    end
    
    xl = XLabelConfig(pc.axesConfig(row,col), label, xGap, positionLabel, 'top', 0); 
    xl.plot();
    
    if i <= 2 && length(CountFig) <= 2
        xl = XLabelConfig(pc.axesConfig(row,col), XLabel, xGap, 'middle', 'bottom', 0);  
        xl.plot
    end
    
    if row == nRows 
        xl = XLabelConfig(pc.axesConfig(row,col), XLabel, xGap, 'middle', 'bottom', 0);  
        xl.plot
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(stPlot,'slev')
        xPositions  = pc.axesConfig(row,col).ax_xLim; 
        yPositions  = [stPlot.slev stPlot.slev];
        axesConfig  = pc.axesConfig(row,col);
        l1          = LineCurve(axesConfig, xPositions, yPositions, '--', [0.75 0.75 0.75], 1, '');
        l1.plot();
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ac = pc.axesConfig(row,col);
    if row > 1 || col > 1
        tl = TextLabel(ac, ID, [0.1 0.1], 'normalized', 0, '', fntsz);
        tl.plot();
    end
    
    h = gcf;
    set(gcf, 'PaperPositionMode','auto')
    
end