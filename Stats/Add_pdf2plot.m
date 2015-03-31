function [PDF yi] = Add_pdf2plot(x, Handle)
% function [PDF yi] = Add_pdf2plot(x, Handle)
%
% 1. Description:
%   Add probability density function on the right of a plot. The variable
%   'x' has to be plotted in the figure identified by 'Handle'. If you do
%   not specify a handle, then this function will assume that x is plotted 
%   in the most recent figure.
% 
%   TO DO: Put mean and deviation. See DSP Guide, Ch 2, page 14
% 
%   Handle - Figure handle
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   N = 20;
%   x = wgn(1000,1,1);
%   Add_pdf2plot(x);
% 
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on  : 21/05/2014
% Last update : 21/05/2014 % Update this date manually
% Last used on: 28/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    figure;
    plot(x);
    grid on
    xlabel('[samples]')
    ylabel('Amplitude')
    Handle      = gcf;
end
Handle_axis = gca;

N = 20; % PDF-resolution
[PDF yi] = Probability_density_function(x,N);

delta_x = max(get(Handle_axis,'XLim'))*0.1;
Xlimits = get(Handle_axis,'XLim');

off_setx = Xlimits(2); 
set(Handle_axis, 'XLim',[Xlimits(1) Xlimits(2)+delta_x])

figure(Handle);
hold on
plot(delta_x*PDF/max(PDF)+off_setx,yi,'r')
legend('Series','pdf')

% rectangle('Position',[-50 0.01 50 20-0.05],'FaceColor','w'), hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
