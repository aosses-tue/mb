function Saveas(h, filename, option)
% function Saveas(h, filename, option)
%
% 1. Description:
%       Function to save figures in vectorial formats. EPS is more suitable 
%       for LaTeX documents and EMF is more suitable if figures are going 
%       to be included in a microsoft PPT presentation.
%           - option.format = 'epsc'; % for eps in color
%           - option.format = 'emf'; for emf format
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium 2014
% Created in    : January-March 2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 20/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    option = [];
end
option = Ensure_field(option,'bPrint',1);
option = ef(option,'format','epsc');

set(h,'PaperType', 'A4')
Check_figure_size(h); % Constrains figure width to A4 Paper size
set(h,'PaperPositionMode', 'auto')
if option.bPrint
%     figure(h);
%     Print_date_on_figure;
end

saveas(h, filename, option.format);

disp(['Figure saved as: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
