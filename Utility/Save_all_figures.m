function Save_all_figures(Handles, outputdir, counter_i)
% function Save_all_figures(Handles, outputdir, counter_i)
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
% Last update on: 26/11/2014 % Update this date manually
% Last use on   : 20/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    Handles = findobj('Type','figure');
end

if nargin < 2
    outputdir = Get_TUe_paths('outputs');
end

if nargin < 3
    counter_i = 1;
end

for i = 1:length(Handles)
    Saveas(Handles(i), [outputdir 'fig-' num2str(i+counter_i-1)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
