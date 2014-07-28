function handle_printed_text = Print_text_on_plot(text2print,tstart,ystart)
% function handle_printed_text = Print_text_on_plot(text2print,tstart,ystart)
%
% 1. Description:
%       Print text on plot, is intended to be a kind of watermark, to label
%       plots in case they are generated automatically 
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 4/6/2014
% Last update: 4/6/2014 % Update this date manually
% Last used: 4/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    ha = gca;
    YLim = get(ha,'YLim');
    ystart = YLim(1) + 0.1*(YLim(2)-YLim(1));
end

if nargin < 2
    ha = gca;
    XLim = get(ha,'XLim');
    tstart = XLim(2);
    XLim(2) = XLim(2) + 0.1*(XLim(2)-XLim(1)); % 10% Extension 
    set(ha,'XLim',XLim)
end
% ystart = 0.2

Color = [0.75 0.75 0.75];
handle_printed_text = text(tstart,ystart, text2print, 'HorizontalAlignment','left','Color',Color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end