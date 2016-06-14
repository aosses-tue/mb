function [filename outsig] = AM_sine_savewave(fc,dur,fs,fmod,Mdept,SPL,dir)
% function [filename outsig] = AM_sine_savewave(fc,dur,fs,Fmod,Mdept,SPL,dir)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fc = 1000;
%       dur = 4;
%       fs = 44100;
%       AM_sine(fc,dur,fs);
%
%       fc  = 1000;
%       dur = 4;
%       fs  = 44100;
%       fmod = 4; % Hz
%       Mdept = 1;
%       SPL = 70;
%       y   = AM_sine(fc,dur,fs,fmod,Mdept,SPL);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 02/04/2015
% Last update on: 02/04/2015 
% Last use on   : 20/05/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7
    try
        dir = Get_TUe_paths('outputs'); % Alejandro's directory
    catch
        dir = pwd;
    end
end

dBFS = 100;

if nargin < 6
    SPL = 65;
end

if nargin < 5
    Mdept = 100;
end

if nargin < 4
    fmod = 4;
end

if nargin < 3
    fs = 44100; % sampling frequency
end

if nargin < 2
    Tmod = 3/fmod;
    dur = 1*Tmod; 
end
        
sig = Create_sin(fc,dur,fs,0);

option = 'm';
start_phase = -pi/2; % Use -pi/2 to begin in minimum
outsig = ch_am(sig,fmod,Mdept,option,fs,start_phase);

outsig = setdbspl(outsig,SPL,dBFS); % same than applying calibration factor: cal = From_dB(-dBFS)*( From_dB(SPL)/mean(rms(y)) );

% sound(y,fs);
filename = [dir sprintf('AM-tone-fc-%.0f_fmod-%.0f_mdept-%.0f-SPL-%.0f-dB',fc,fmod,Mdept,SPL)];
if nargout < 2
    Wavwrite(outsig,fs,filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
