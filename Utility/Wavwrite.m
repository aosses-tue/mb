function Wavwrite(x,fs,outputfilename)
% function Wavwrite(x,fs,outputfilename)
% 
% 1. Description:
%       Write Microsoft WAVE (".wav") sound file
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 4. Stand-alone example:
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on    : 13/05/2014
% Last update on: 13/05/2014 % Update this date manually
% Last use on   : 16/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavwrite(x,fs,[outputfilename]);

disp([mfilename '.m: file ' outputfilename ' created'])
try
    fprintf('%s.m: RMS of %.2f [dBFS]\n',mfilename,rmsdb(x))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
