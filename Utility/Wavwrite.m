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
% Created on: 13/5/2014
% Last update: 13/5/2014 % Update this date manually
% Last used: 27/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavwrite(x,fs,[outputfilename]);

disp([mfilename '.m: file ' outputfilename ' created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end