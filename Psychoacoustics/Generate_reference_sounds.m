function y = Generate_reference_sounds
% function y = Generate_reference_sounds
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc  = 1000;
T   = 1/fc;
fmod= 70;
Tmod= 3/fmod;
dur = 1;
fs  = 44100;
sig = Create_sin(fc,dur,fs,0);
m   = 100;
option = 'm';
start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
ch_am(sig,fmod,m,option,fs,start_phase);
y = ch_am(sig,fmod,m,option,fs,start_phase);

lvl     = 70;
lvlAMT  = lvl + 10;
y = setdbspl(y,lvlAMT);

Wavwrite(y,fs,[Get_TUe_paths('outputs') 'ref_sharpness']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
