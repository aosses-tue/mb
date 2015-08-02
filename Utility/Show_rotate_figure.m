function Show_rotate_figure(delay_time,ha,h)
% function Show_rotate_figure(delay_time,ha,h)
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
% Created on    : 25/07/2015
% Last update on: 25/07/2015 
% Last used on  : 25/07/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    delay_time = 1; % 1 second
end

if nargin < 1
    ha = gca;
end

if nargin < 2
    h = gcf;
end

positions = -50:10:0;
for i = 1:length(positions)
    set(ha,'View',[positions(i) 50])
    figure(h);
    title(sprintf('Az, el: %.1f, %.1f',positions(i),50))
    pause(delay_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])