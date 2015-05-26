function [numH7, denH7, numH14, denH14, numH30, denH30, numH36, denH36, numH66, denH66] = fir_hweight(Fs)
% function [numH7, denH7, numH14, denH14, numH30, denH30, numH36, denH36, numH66, denH66] = fir_hweight(Fs)
%
% 1. Description:
%       This file calculates the coefficients of the FIR filters that are
%       an approximation of the weighting function Hi(fmod) used by Daniel
%       and Weber
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Adapted by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Author        : Vincent Jourdes
% Created on    : 26/05/2015
% Last update on: 26/05/2015 % Update this date manually
% Last use on   : 26/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    Fs = 40960;
end

% Degree of the filter
nn = 4;
dd = 6;

% Modification of the degree for signals with a sampling frequency of 40960 Hz
if Fs == 40960
    nn = 3;
    dd = 6;
else
    nn = 4;
    dd = 6;
end

% H7
f7 = [0 17 23 25 32 37 48 67 90 114 171 206 247 294 358 500 Fs/2]/(Fs/2);
H7 = [0 0.8 0.95 0.975 1 0.975 0.9 0.8 0.7 0.6 0.4 0.3 0.2 0.1 0 0 0];
edges7 = [0 32 37 Fs/2]/(Fs/2);
w = [1 1 1 3 7 3 1 1 1 1 1 1 1 1 1 1 1];
[numH7, denH7] = iirlpnorm(nn,dd,f7,edges7,H7,w);   % calculate the coefficients
                                                    % of a fulter that is an 
                                                    % approximation of H7

% H14
f14 = [0 32 43 56 69 92 120 142 165 231 277 331 397 502 1000 Fs/2]/(Fs/2);
H14 = [0 0.8 0.95 1 0.975 0.9 0.8 0.7 0.6 0.4 0.3 0.2 0.1 0 0 0];
edges14 = [0 56 69 Fs/2]/(Fs/2);
w = [1 1 2 8 2 1 1 1 1 1 1 1 1 1 1 1];
[numH14, denH14] = iirlpnorm(nn,dd,f14,edges14,H14,w);

% H30
if Fs == 44100
    nn = 3;
    dd = 6;
else
    nn = 4;
    dd = 6;
end

f30 = [ 0 23.5 34 47 56 63 72 78 87 100 115 135 159 172 194 215 244 290 ...
        348 415 500 645 1000 Fs/2]/(Fs/2);
H30 = [0 0.4 0.6 0.8 0.9 0.95 0.98 1 .99 .975 .95 .9 .85 .8 .7 .6 .5 .4 .3 .2 .1 0 0 0];
edges30 = [0 78 87 Fs/2] / (Fs/2);
w = [1 1 1 1 1 1 4 10 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
[numH30, denH30] = iirlpnorm(nn,dd,f30,edges30,H30,w);

% H36
f36 = [ 0 19 44 52.5 58 75 80 90 101.5 114.5 132.5 143.5 165.5 197.5 241 290 ...
        348 415 500 645 1000 Fs/2]/(Fs/2);
H36 = [0 0.4 0.8 0.9 0.95 1 .98 .97 .95 .9 .85 .8 .7 .6 .5 .4 .3 .2 .1 0 0 0];
edges36 = [0 75 80 Fs/2]/(Fs/2);
w = [1 1 1 1 1 6 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
[numH36, denH36] = iirlpnorm(nn,dd,f36,edges36,H36,w);

% H66
f66 = [ 0 15 41 49 53 64 71 88 94 106 115 137 180 238 290 348 415 500 645 ... 
        1000 Fs/2]/(Fs/2);
H66 = [0 0.4 0.8 0.9 0.965 0.99 1 0.95 0.9 0.85 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 0 0];
edges66 = [0 71 88 Fs/2]/(Fs/2);
w = [1 1 1 1 1 4 8 4 1 1 1 1 1 1 1 1 1 1 1 1 1];
[numH66, denH66] = iirlpnorm(nn,dd,f66,edges66,H66,w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
