function Save_figure_as(h, filename, format)
% function Save_figure_as(h, filename, format)
%
% 1. Description:
%       Function to save figures in vectorial formats. EPS is more suitable 
%       for LaTeX documents and EMF is more suitable if figures are going 
%       to be included in a microsoft PPT presentation.
%           - format = 'epsc'; % for eps in color
%           - format = 'emf'; for emf format
%       Other formats might be available
%
% 2. Stand-alone example:
%       dur = 10e-3;
%       fs = 44100;
%       y  = wgn(round(dur*fs),1,1);
%       t = ( 1:length(y) )/fs;
%       figure;
%       plot(t*1000,y)
%       xlabel('time [ms]'), ylabel('amplitude'), title('white noise')
%       Save_figure_as(gcf,'white-noise','epsc') % to save it in eps color
%       Save_figure_as(gcf,'white-noise','emf') % to save it in emf
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium 2014
% Created in    : January-March 2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 28/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    format = 'epsc';
end

str = fileparts(filename);

if strcmp(str,'') % if not directory is specified...
    str = [pwd delim];
    filename = [str filename];
end

set(h,'PaperType', 'A4')
set(h,'PaperPositionMode', 'auto')

saveas(h, filename, format);

disp(['Figure saved as: ' filename ', in ' format ' format'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
