function [y t] = Create_sin(f,dur,fs,window)
% function [y t] = Create_sin(f,dur,fs,window)
%
% 1. Description:
%   f = frequency [Hz], default = 1000 Hz
%   dur = duration [s], default = 1 s
%   fs = sampling frequency, default = 44100 Hz
%   window = Window to be applied to signal
%       0 = none
%       1 = hanning
%       2 = triangular
%
% 2. Stand-alone example, sine wave with default values:
%   f       = 1000;
%   dur     = 20e-3;
%   [y, t]  = Create_sin(f,dur);
%   plot(t,y); xlabel('time [s]'); ylabel('Amplitude')
% 
% 3. Additional info:
%   Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 26/05/2014
% Last update on: 26/06/2014 
% Last use on   : 16/03/2016 
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
    window = 0;
end

t = 0: 1/fs : dur - 1/fs;
t = t';

w = 2*pi*f;
y =  sin(w*t);
T = 1/f;

try
    [win, wtype] = Get_window(window,y);
    y = y.*win;
catch
    if window ~= 0
        error('The window has not been applied, use window = 0 (rectangular) otherwise get the script Get_window');
    end
end

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