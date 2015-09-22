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
% Last use on   : 28/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    option = [];
else
    if ischar(option)
        option.format = option;
    else
        option = ef(option,'format','epsc');
    end
end

str = fileparts(filename);

if strcmp(str,'')
    try
        str = Get_TUe_paths('outputs');
    catch
        str = [pwd delim];
    end
    filename = [str filename];
end

option = Ensure_field(option,'bPrint',1);
option = Ensure_field(option,'bScale',0);


set(h,'PaperType', 'A4')

if option.bScale == 1
    try
        Check_figure_size(h); % Constrains figure width to A4 Paper size
    end
end

set(h,'PaperPositionMode', 'auto')
if option.bPrint
%     figure(h);
%     Print_date_on_figure;
end

saveas(h, filename, option.format);

disp(['Figure saved as: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
