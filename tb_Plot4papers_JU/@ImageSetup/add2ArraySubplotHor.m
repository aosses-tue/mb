function obj = add2ArraySubplotHor(obj);
% function obj = add2ArraySubplotHor(obj);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj.figs2subplotsArray(obj,'Handles',obj.arrayAddedHandles,...
                          'Tiling',obj.I_Matrix,...
                          'Direction','horizontal');
obj.prepareFigures('FigHandle',obj.hOutFig);