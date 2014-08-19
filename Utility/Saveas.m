function Saveas(h, filename, option)
% function Saveas(h, filename, option)
%
% Function to save EPS files
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, Belgium 2014
% Created in    : January-March 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 18/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    option = [];
end
option = Ensure_field(option,'bPrint',1);

set(h,'PaperType', 'A4')
Check_figure_size(h); % Constrains figure width to A4 Paper size
set(h,'PaperPositionMode', 'auto')
if option.bPrint
%     figure(h);
%     Print_date_on_figure;
end

saveas(h, filename, 'epsc');

disp(['Figure saved as: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end