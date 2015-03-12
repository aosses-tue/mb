function hOutput = Figure2paperfigureT2(nHandle,numRows,numCols,plotOpts)
% function hOutput = Figure2paperfigureT2(nHandle,numRows,numCols,plotOpts)
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 26/09/2014
% Last update on: 28/02/2015 % Update this date manually
% Last use on   : 28/02/2015 % Update this date manually
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
plotOpts = ef(plotOpts,'I_Width'   , 11);
plotOpts = ef(plotOpts,'I_Height'  , 10);

I_FontSize = plotOpts.I_FontSize;
I_Width    = plotOpts.I_Width;
I_Height   = plotOpts.I_Height;

hM = ImageSetup; 
hM.I_Matrix      = [numRows,numCols];
hM.I_FontSize    = I_FontSize; 
hM.I_FontName    = 'Arial'; 
hM.I_Width       = I_Width;
hM.I_Height      = I_Height;
hM.I_TitleInAxis = 1;
hM.I_Space       = [0.05,0.05];
hM.I_Legend      = 'off';

% hM.I_Ylim       = [-50,0]; % Uncomment for fixing the limits in the y-axis
% hM.I_Xlim       = [0 24];

hM.I_Grid        = 'on'; 
hM.I_KeepColor   = 1; 
hM.I_Handles     = nHandle;
hM.prepareAllFigures;
hM.arrayAddedHandles = 1;

add2ArraySubplotVer(hM);

hOutput = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end