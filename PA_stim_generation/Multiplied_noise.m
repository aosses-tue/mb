function y = Multiplied_noise(fc,BW,SPL,dur,fs)
% function y = Multiplied_noise(fc,BW,SPL,dur,fs)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fc = 1300;
%       BW = 100;
%       SPL = 60;
%       dur = 1;
%       fs = 44100;
%       Multiplied_noise(fc,BW,SPL,dur,fs);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/08/2015
% Last update on: 20/08/2015 
% Last use on   : 20/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = Create_sin(fc,dur,fs);
n = AM_random_noise(0,BW/2,60,dur,fs,0);

y = s.*n;
y = setdbspl(y,SPL); % calibration

if nargout == 0
    
    filename = [Get_TUe_paths('outputs') sprintf('MN-fc-%.0f_BW-%.0f_SPL-%.0f-dB',fc,BW,SPL)];
    Wavwrite(y,fs,filename);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
