function hOutput = Figure2paperfigureT(nHandle,numRows,numCols,plotOpts)
% function hOutput = Figure2paperfigureT(nHandle,numRows,numCols,plotOpts)
%
% 1. Description:
%       Adapts figure in nHandle to be paper-compatible, using Jaime Undurraga's
%       GUI. hOutput is the output figure handle of the resulting figure
% 
% 2. Stand-alone example:
%  
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 26/09/2014
% Last update on: 24/07/2015 
% Last use on   : 24/07/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    numRows = 4;
end

if nargin < 3
    numCols = 2;
end

if nargin < 4
    plotOpts = [];
end

plotOpts = ef(plotOpts,'I_FontSize',14);
plotOpts = ef(plotOpts,'I_Width'   ,10);
plotOpts = ef(plotOpts,'I_Height'  , 8);
plotOpts = ef(plotOpts,'I_TitleInAxis',1);
plotOpts = ef(plotOpts,'I_KeepColor',1);
plotOpts = ef(plotOpts,'I_KeepTicks',1);
plotOpts = ef(plotOpts,'I_Legend','off');
plotOpts = ef(plotOpts,'bAddVertical',1);

bAddVertical = plotOpts.bAddVertical;

I_FontSize = plotOpts.I_FontSize;
I_Width    = plotOpts.I_Width;
I_Height   = plotOpts.I_Height;
I_TitleInAxis = plotOpts.I_TitleInAxis;
I_KeepColor = plotOpts.I_KeepColor;
I_KeepTicks = plotOpts.I_KeepTicks;
I_Legend    = plotOpts.I_Legend;

hM = ImageSetup; 
hM.I_Matrix      = [numRows,numCols];
hM.I_FontSize    = I_FontSize; 
hM.I_FontName    = 'Arial'; 
hM.I_Width       = I_Width;
hM.I_Height      = I_Height;
hM.I_TitleInAxis = I_TitleInAxis;
hM.I_Space       = [0.05,0.05];
hM.I_Legend      = I_Legend;

% hM.I_Ylim       = [-50,0]; % Uncomment for fixing the limits in the y-axis
% hM.I_Xlim       = [0 24];

if isfield(plotOpts,'I_Ylim')
    hM.I_Ylim    = plotOpts.I_Ylim; 
end
if isfield(plotOpts,'I_Xlim')
    hM.I_Xlim    = plotOpts.I_Xlim; 
end

hM.I_Grid        = 'on'; 
hM.I_KeepColor   = I_KeepColor; 
hM.I_KeepYTicks  = I_KeepTicks;
hM.I_Handles     = nHandle;
hM.prepareAllFigures;
hM.arrayAddedHandles = 1;

if bAddVertical
    add2ArraySubplotVer(hM);
else
    add2ArraySubplotHor(hM);
end

hOutput = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
