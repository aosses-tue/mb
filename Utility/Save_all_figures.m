function Save_all_figures(Handles)
% function Save_all_figures(Handles)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       figure;
%       figure;
%       figure;
%       Save_all_figures;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/08/2014
% Last update on: 12/08/2014 % Update this date manually
% Last use on   : 12/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    Handles = findobj('Type','figure');
end

for i = 1:length(Handles)
    Saveas(Handles, [Get_TUe_paths('outputs') 'fig-' num2str(i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
