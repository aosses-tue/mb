function [b,a] = hpf(fc,fs)
% function [b,a] = hpf(fc,fs)
%
% 1. Description:
%       Designs a Kaiser window with cut-off at fc, sampled at fs
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 5/6/2014
% Last update: 5/6/2014 % Update this date manually
% Last used: 5/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
% F - band edges
% A - desired amplitudes at band edges
% DEV - is a vector of maximum deviations or ripples
% Fs - sampling freq

if nargin < 1
    fc = 50;
end

if nargin < 2
    fs = 10000;
end
    
DEV = 2*[0.01 0.1];

[n,Wn,beta] = kaiserord( [0.7*fc fc], [0 sqrt(2)], DEV, fs );

b = fir1(n,Wn,'high',kaiser(n+1,beta));

% freqz(b,1,512)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end