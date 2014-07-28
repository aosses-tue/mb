function Saveas(h, filename)
% function Saveas(h, filename)
%
% Function to save EPS files
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, Belgium 2014
% Created in: January-March 2014
% Last update: 24/07/2014 % Update this date manually
% Last used: 24/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h,'PaperType', 'A4')
Check_figure_size(h); % Constrains figure width to A4 Paper size
set(h,'PaperPositionMode', 'auto')

saveas(h, filename, 'epsc');

disp(['Figure saved as: ' filename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end