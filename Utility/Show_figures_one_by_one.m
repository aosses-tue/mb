function figHandles = Show_figures_one_by_one(delay_time)
% function figHandles = Show_figures_one_by_one(delay_time)
%
% 1. Description:
%       Show figures at a time span of delay_time seconds
% 
% 2. Stand-alone example:
%       % Make sure you have several figures open...
%       Show_figures_one_by_one(1);
% 
% 3. Additional info:
%   Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/06/2014
% Last update on: 05/06/2014 % Update this date manually
% Last used on  : 25/27/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    delay_time = 1; % 1 second
end

figHandles = findobj('Type','figure');

for i = length(figHandles):-1:1
    figure(figHandles(i));
    pause(delay_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])