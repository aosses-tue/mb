function [b,a] = lpf(fc,fs)
% function [b,a] = lpf(fc,fs)
%
% 1. Description:
%       Designs a Kaiser window (LPF) with cut-off at fc, sampled at fs
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

if nargin < 2
    fs = 10000;
end

if nargin < 1
    fc = fs/2-1;
end

if fc >= fs/2
    fc = fs/2-1;
    f2 = fc;
    f1 = 0.9*fc;
else
    f2 = min(1.1*fc,fs/2 - 1);
    f1 = 0.9*f2; 
end

DEV = 2*[0.01 0.1];

[n,Wn,beta] = kaiserord( [f1 f2], [sqrt(2) 0], DEV, fs );

b = fir1(n,Wn,'low',kaiser(n+1,beta));

% freqz(b,1,512)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end