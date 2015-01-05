function y = Do_zeropadding_in_samples(x,si,st)
% function y = Do_zeropadding_in_samples(x,si,st)
%
% 1. Description:
%       x - signal to be zero padded
%       si - intial sample
%       st - total samples (length of the output signal)
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/01/2015
% Last update on: 03/01/2015 % Update this date manually
% Last use on   : 03/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ni = length(x);

y = [zeros(si,1); x; zeros(st-si-Ni,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
