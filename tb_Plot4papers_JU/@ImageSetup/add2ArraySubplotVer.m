function obj = add2ArraySubplotVer(obj);
% function obj = add2ArraySubplotVer(obj);
%
% Programmed by Jaime Undurraga, modified by Alejandro Osses
% Last use on   : 11/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj.figs2subplotsArray(obj,'Handles',obj.arrayAddedHandles,...
    'Tiling',obj.I_Matrix,'Direction','vertical');
obj.prepareFigures('FigHandle',obj.hOutFig);