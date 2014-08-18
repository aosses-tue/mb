function [sig t] = fm(fc, dur, fs, fmod, deltaf)
% function [sig t] = fm(fc, dur, fs, fmod, deltaf)
% 
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/08/2014
% Last update on: 16/08/2014 % Update this date manually
% Last use on   : 16/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    fc = 1000;
end

if nargin < 2
    dur = 1;
end

if nargin < 3
    fs = 44100;
end

if nargin < 4
    fmod = 70; 
end

if nargin < 5
    deltaf = 700; 
end

t = 0: 1/fs : dur - 1/fs;
t = t';

w = 2*pi*fc;
start_phase = -pi/2;
sig =  sin(w*t - (deltaf/fmod)*cos(2*pi*fmod*t+start_phase));
% T = 1/fc;

if nargout == 0
    figure;
    plot( t,sig )
    xlabel('time [s]')
    grid on
    xlim([0 2/fmod])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
