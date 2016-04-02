function output = Get_colour(str_colour)
% function output = Get_colour(str_colour)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 21/03/2016
% Last update on: 21/03/2016 
% Last use on   : 21/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch str_colour
    case 'blue'
        output = [0 0 1];
    case 'blue_dark'
        output = [0 0 0.5];
    case 'red'
        output = [1 0 0];
    case 'red_dark'
        output = [0.5 0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
