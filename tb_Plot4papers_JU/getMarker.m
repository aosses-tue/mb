function [marker] = getMarker(i)
% function [marker] = getMarker(i)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Jaime Undurraga
% Created in   : 2008-2012
% Downloaded on: 30/09/2012
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update  : 24/07/2014 % Update this date manually
% Last used    : 24/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Markers = ['s' 'd' '+' 'o' '*' '.' 'x'  '^' 'v' '>' '<' 'p' 'h'];
LM      = length(Markers);
aux     = mod(i-1,LM)+1;
marker  = Markers(aux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end