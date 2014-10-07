function f = bark2hz(z)
% function f = bark2hz(z)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       z = 5; % Bark
%       f = bark2hz(z);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/09/2014
% Last update on: 25/09/2014 % Update this date manually
% Last use on   : 25/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = 1000;
k = -16:12;
f = f0*2.^(k/3); 

zt = hz2bark(f);

f = interp1(zt,f,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
