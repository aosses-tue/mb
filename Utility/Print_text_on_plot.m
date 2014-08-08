function handle_printed_text = Print_text_on_plot(text2print,tstart,ystart,color_scale)
% function handle_printed_text = Print_text_on_plot(text2print,tstart,ystart,color_scale)
%
% 1. Description:
%       Print text on plot, is intended to be a kind of watermark, to label
%       plots in case they are generated automatically 
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       figure;
%       plot([0 0], [1 1])
%       color_scale = 0.25; % gray almost black
%       tstart = 0.2;
%       ystart = 0;
%       text2print = 'Water mark';
%       Print_text_on_plot(text2print,tstart,xstart,color_scale);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 04/06/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    ha = gca;
    YLim = get(ha,'YLim');
    ystart = YLim(1) + 0.05*(YLim(2)-YLim(1));
end

if nargin < 2
    ha = gca;
    XLim = get(ha,'XLim');
    tstart = XLim(1) + 0.3*(XLim(2)-XLim(1));
end

if nargin < 4
    color_scale = 0.75;
end
% ystart = 0.2

Color = color_scale*[1 1 1];
handle_printed_text = text(tstart,ystart, text2print, 'HorizontalAlignment','left','Color',Color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end