function X2 = subtract_dB(X1,delta_dB)
% function X2 = subtract_dB(X1,delta_dB)
%
% 1. Description:
%
% 2. Stand-alone example:
%       X1 = 60;
%       X2 = subtract_dB(X1,1);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/05/2015
% Last update on: 06/05/2015 % Update this date manually
% Last use on   : 06/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if delta_dB == 0
    delta_dB = eps;
end

Y   = X1 + delta_dB;
X2  = To_dB( From_dB(Y)-From_dB(X1)+eps );

if X2 < 0
    X2 = 0.01;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
