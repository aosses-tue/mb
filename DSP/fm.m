function [sig t] = fm(fc, dur, fs, fmod, deltaf)
% function [sig t] = fm(fc, dur, fs, fmod, deltaf)
% 
% 1. Description:
%       It generates an FM tone with centre frequency fc, duration dur [s]
%       at a sample frequency fs. The modulation frequency is fmod and the 
%       frequency deviation is given by +/- deltaf so the frequency of the 
%       generated tone oscilates between fc - deltaf and fc + deltaf.
% 
% 3. Stand-alone example:
%       fc = 1000; % Hz
%       dur = 100e-3; % 100 ms
%       fs = 44100;
%       fmod = 70; % Hz
%       deltaf = 700; % Hz
%       fm(fc,dur,fs,fmod,deltaf);
% 
% 2. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 16/08/2014
% Last update on: 27/11/2014 
% Last use on   : 27/11/2014 
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
if fmod ~= 0
    sig =  sin( w*t - (deltaf/fmod)*cos(2*pi*fmod*t+start_phase) );
else
    sig =  sin( w*t ); % no modulation
end

if nargout == 0
    figure;
    plot( t,sig )
    xlabel('time [s]')
    ylabel('Amplitude')
    title(sprintf('2 modulation periods for FM tone, fc=%.0f [Hz], fmod=%.0f [Hz], fdev=+/-%.2f [Hz]',fc,fmod,deltaf))
    grid on
    xlim([0 2/fmod])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
