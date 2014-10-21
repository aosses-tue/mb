function [y t] = Create_sin_phase(f,start_phase,dur,fs,window)
% function [y t] = Create_sin_phase(f,start_phase,dur,fs,window)
%
% 1. Description:
%   f = frequency [Hz], default = 1000 Hz
%   start_phase = start phase [rad]
%   dur = duration [s], default = 1 s
%   fs = sampling frequency, default = 44100 Hz
%   window = Window to be applied to signal
%       0 = none
%       1 = hanning
%       2 = triangular
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example, sine wave with default values:
%   f       = 1000;
%   start_phase = pi/2; % starts in amplitude 1
%   dur     = 20e-3;
%   [y, t]  = Create_sin_phase(f,start_phase, dur);
%   figure; plot(t,y); xlabel('time [s]'); ylabel('Amplitude')
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/10/2014
% Last update on: 21/10/2014 % Update this date manually
% Last use on   : 21/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    f = 1000;
end

if nargin < 3
    start_phase = 0;
end

if nargin < 3
    dur = 1;
end

if nargin < 4
    fs = 44100;
end

if nargin < 5
    window = 0;
end

t = 0: 1/fs : dur - 1/fs;
t = t';

w = 2*pi*f;
y = sin(w*t + start_phase);
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