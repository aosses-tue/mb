function y = dB2Pa(x)
% function y = dB2Pa(x)
%
% 1. Description:
%       Converts from Sound pressure level x [dB] into pressure y [Pa].
% 
% 2. Stand-alone example:
%       y = dB2Pa(60);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/05/2015
% Last update on: 06/05/2015 % Update this date manually
% Last use on   : 06/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0 = 2e-5;
y = p0 * From_dB(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
