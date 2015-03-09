function d = m2d(m)
% function d = m2d(m)
%
% 1. Description:
%       Converts modulation factor 'm' to modulation depth 'd'. Enter a 
%       value 'm' between 1 and 100 [%]
% 
% 2. Stand-alone example:
%       d = m2d(94)
% 
% 3. Additional info:
%       Tested cross-platform: yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/08/2014
% Last update on: 16/08/2014 % Update this date manually
% Last use on   : 05/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if m > 1
    m = m/100;
end

d = 20 * log10( (1+m)/(1-m) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
