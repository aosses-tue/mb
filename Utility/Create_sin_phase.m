function [y t] = Create_sin_phase(f,dur,fs,start_phase)
% function [y t] = Create_sin_phase(f,dur,fs,start_phase)
%
% 1. Description:
%   Sine with initial phase, no window will be applied.
%       f = frequency [Hz], default = 1000 Hz
%       dur = duration [s], default = 1 s
%       fs = sampling frequency, default = 44100 Hz
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example, sine wave with default values:
%   f       = 1000;
%   dur     = 1;
%   start_phase = 0;
%   fs      = 44100;
%   [y, t]  = Create_sin_phase(f,dur,fs,start_phase);
%   plot(t,y); 
%   xlabel('time [s]'); 
%   ylabel('Amplitude')
%   m       = 100;
%   fmod    = 70;
%   option  = 'm';
%   start_phase = pi/2; % modulation starts in maximum
%   ch_am(y,fmod,m,option,fs,start_phase); 
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 26/05/2014
% Last update on: 26/06/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    f = 1000;
end

if nargin < 2
    dur = 1;
end

if nargin < 3
    fs = 44100;
end

if nargin < 4
    start_phase = pi/2; % max
end

window = 0;

t = 0: 1/fs : dur - 1/fs;
t = t';

w = 2*pi*f;
y =  sin(w*t+start_phase);
T = 1/f;

[win, wtype] = Get_window(window,y);

y = y.*win;

if nargout == 0
    plot(t,y);
    xlim([0 2*T])
    title(['Sine wave, 2 periods, ' wtype ' window'])
    xlabel('time [s]')
    ylabel('Amplitude')
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end