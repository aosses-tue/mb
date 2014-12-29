function obj = prepareAllFigures(obj)
% function obj = prepareAllFigures(obj)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% replot subplots figures
    if obj.updateContainerFigure
        if ishandle(obj.hOutFig)
            delete(obj.hOutFig)
        end
        if ~isempty(obj.verAddedHandles)
            obj.addVerHandle([]);
        end
        if ~isempty(obj.horAddedHandles)
            obj.addHorHandle([]);
        end
        if ~isempty(obj.arrayAddedHandles)
            obj.addHandle2Array([]);
        end
        obj.updateContainerFigure = false;
    end
    %%
    if ~isprop(obj,'I_Handles')
        obj.findImageHandles;
    end
    
    for i=1:length(obj.I_Handles)
        obj.prepareFigures('FigHandle', obj.I_Handles(i));
    end
end