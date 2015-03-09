function m = d2m(d)
% function m = d2m(d)
%
% 1. Description:
%       Converts modulation depth 'd' into modulation factor 'm'. Enter a 
%       value 'd' in [dB]
% 
% 2. Stand-alone example:
%       m = d2m(40)
% 
% 3. Additional info:
%       Tested cross-platform: yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/03/2015
% Last update on: 05/03/2015 % Update this date manually
% Last use on   : 05/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = ( 10.^(d/20)-1 )/( 10.^(d/20)+1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
