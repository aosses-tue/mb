function y = AM_sine(fc,dur,fs,fmod,Mdept,SPL,dBFS)
% function y = AM_sine((fc,dur,fs,Fmod,Mdept,SPL,dBFS)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fc = 1000;
%       dur = 4;
%       fs = 44100;
%       AM_sine(fc,dur,fs);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 02/04/2015
% Last update on: 02/04/2015 % Update this date manually
% Last use on   : 02/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin < 7
%     dBFS = 100;
% end

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
start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
y = ch_am(sig,fmod,Mdept,option,fs,start_phase);
    
% RMS = 65-100;
y = setdbspl(y,SPL);

if nargout == 0
    
    sound(y,fs);
    filename = [Get_TUe_paths('outputs') sprintf('tone-fc-%.0f__Fmod-%.0f_Mdept-%.0f',fc,fmod,Mdept)];
    Wavwrite(y,fs,filename);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
