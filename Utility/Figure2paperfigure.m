function hOutput = Figure2paperfigure(nHandle,numRows,numCols)
% function hOutput = Figure2paperfigure(nHandle,numRows,numCols)
%
% 1. Description:
%       Adapts figure in nHandle to be paper-compatible, using Jaime Undurraga's
%       GUI. hOutput is the output figure handle of the resulting figure
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 26/09/2014
% Last update on: 26/09/2014 % Update this date manually
% Last use on   : 07/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    numRows = 4;
end

if nargin < 3
    numCols = 2;
end

hM = ImageSetup; 
hM.I_Matrix      = [numRows,numCols];
hM.I_FontSize    = 12; 
hM.I_FontName    = 'Arial'; 
hM.I_Width       = 8;
hM.I_Height      = 8;
hM.I_TitleInAxis = 1;
hM.I_Space       = [0.01,0.01];

% hM.I_Ylim       = [-50,0]; % Uncomment for fixing the limits in the y-axis
% hM.I_Xlim       = [0 24];

hM.I_Grid = 'on'; 
hM.I_KeepColor = 1; 
hM.I_Handles = nHandle;
hM.prepareAllFigures;
hM.arrayAddedHandles = 1;

add2ArraySubplotVer(hM);

hOutput = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
