function [w,t] = Conv(u,v,info)
% function [w,t] = Conv(u,v,info)
%
% 1. Summary:
%   Convolution and polynomial multiplication. Same as conv, but gives 
%   possibility to get time.
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Use:
%   Explanation of this function
%
% 4. Stand-alone example:
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 12/5/2014
% Last update: 12/5/2014 % Update this date manually
% Last used: 12/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


info    = Ensure_field(info,'fs',48000);

t       = ( 0:length(u)+length(v)-1 )/info.fs;
t       = t(:); % Forces column vector

w       = conv(u, v); % M = N+N-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])