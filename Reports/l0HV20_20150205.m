function y = l0HV20_20150205
% function y = l0HV20_20150205
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/02/2015
% Last update on: 05/02/2015 % Update this date manually
% Last use on   : 05/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 4; % input('geef aantal segmenten: ');
x_coor = [0 0]; %input('geed x coordinaat oorsprong: ');
y_coor = [0 0]; %input('geed y coordinaat oorsprong: ');
r = 2;%input('geef de straal: ')
plotPolygon(n,x_coor,y_coor,r)

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPolygon(n,x,y,r)

theta = 2*pi/n;
m = [cos(theta) sin(theta); -sin(theta) cos(theta)];
xx(1) = r;
yy(1) = 0;

for i = 2:n+1
    u = m*[xx(i-1) yy(i-1)]';
    xx(i) = u(1);
    yy(i) = u(2);
end

plot(x + xx, y+yy)
axis([x-r x+r y-r y+r])
