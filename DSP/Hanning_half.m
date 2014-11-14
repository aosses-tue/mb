function w = Hanning_half(N)
% function w = Hanning_half(N)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/11/2014
% Last update on: 13/11/2014 % Update this date manually
% Last use on   : 13/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    N = 8192;
    warning('Using default window length...')
end

Hann_tmp = hanning(N, 'periodic'); 
w = Hann_tmp;
w(N/2+1:end) = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
