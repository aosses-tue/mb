function obj = add2ArraySubplotVer(obj);
% function obj = add2ArraySubplotVer(obj);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj.figs2subplotsArray(obj,'Handles',obj.arrayAddedHandles,...
    'Tiling',obj.I_Matrix,'Direction','vertical');
obj.prepareFigures('FigHandle',obj.hOutFig);