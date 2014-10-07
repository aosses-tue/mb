function z = hz2bark(f)
% function z = hz2bark(f)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       f = 1000;
%       z = hz2bark(f)
%
% 4. Reference:
%       Osses2010
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 15/08/2014
% Last update on: 15/08/2014 % Update this date manually
% Last use on   : 18/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 13*atan( 0.76*(f/1000) ) + 3.5*atan( (f/(1000*7.5)).^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
