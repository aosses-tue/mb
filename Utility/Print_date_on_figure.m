function y = Print_date_on_figure
% function y = Print_date_on_figure
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       figure;
%       plot([0 0], [1 1])
%       Print_date_on_figure;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    ha = gca;
catch
    warning('No figure/axis handle detected...')
end

p = Get_date;

txt2print = ['Printed on ' p.date2print];

if isunix
    txt2print = [txt2print ' by AO-u'];
else
    txt2print = [txt2print ' by AO-w7'];
end
    
Print_text_on_plot(txt2print);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
