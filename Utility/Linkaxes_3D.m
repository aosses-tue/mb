function Linkaxes_3D(handles)
% function Linkaxes_3D(handles)
%
% 1. Description:
%       Link axes for 3D plots
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 1/07/2014
% Last update on: 1/07/2014 % Update this date manually
% Last used on  : 1/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ishghandle(handles)
    
   hlink    = linkprop(handles,{'CameraPosition','CameraUpVector'});
   key      = 'graphics_linkprop';
   % Store link object on first subplot axes
   setappdata(handles(1),key,hlink);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end