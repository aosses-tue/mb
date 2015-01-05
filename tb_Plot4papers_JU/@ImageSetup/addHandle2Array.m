function obj = addHandle2Array(obj,newFigHandle,Direction)
obj.checkFigHandles;
obj.arrayAddedHandles = [obj.arrayAddedHandles newFigHandle];
obj.arrayAddedHandles = obj.unique1(obj.arrayAddedHandles);
cbar_handle = findobj(newFigHandle,'tag','Colorbar');
if strcmp(Direction,'horizontal')
    obj.add2ArraySubplotHor;
else
     obj.add2ArraySubplotVer;
end
