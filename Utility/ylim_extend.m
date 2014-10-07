function [y_new x_new] = ylim_extend(hAxis,Ratio)
% function [y_new x_new] = ylim_extend(hAxis,Ratio)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 29/09/2014
% Last update on: 29/09/2014 % Update this date manually
% Last use on   : 29/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    Ratio = 1.25;
end

x_curr = get(hAxis,'XLim');
y_curr = get(hAxis,'YLim');
Ymax = Ratio*y_curr(2);

y_new = [y_curr(1) Ymax];
x_new = x_curr;
set(hAxis,'YLim',y_new);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
