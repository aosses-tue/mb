function [w, corr_dB] = Hanning_half(N)
% function [w, corr_dB] = Hanning_half(N)
%
% 1. Description:
%       Half hanning window. The output w is already compensated with a
%       gain of corr_dB.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/11/2014
% Last update on: 20/11/2014 
% Last use on   : 20/11/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    N = 8192;
    warning('Using default window length...')
end

Hann_tmp = hanning(N, 'symmetric'); 
w = Hann_tmp;
w(N/2+1:end) = 0; 

corr_dB = To_dB(1/(length(w)*mean(w)));
w = From_dB(corr_dB)*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
