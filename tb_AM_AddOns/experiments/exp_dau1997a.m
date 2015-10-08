function masker = exp_dau1997a(fs,nFig)
% function masker = exp_dau1997a(fs,nFig)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/10/2015
% Last update on: 06/10/2015 
% Last use on   : 06/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

% f = 1000;
% durramp = 5; % ms
% lvl = 75; % arbitrary

if nargin < 2
    nFig = 3; % modulation detection
end

switch nFig
    case 3
    % To recreate Fig.3

        BW = 3;
        fc = 5000;
        dur = 10;
        SPL = 65;
        masker = AM_random_noise(fc-BW/2,fc+BW/2,SPL+3,dur,fs); 
    
    case 4
    % To recreate Fig.4

        BW = 31;
        fc = 5000;
        dur = 10;
        SPL = 65;
        masker = AM_random_noise(fc-BW/2,fc+BW/2,SPL+3,dur,fs);         
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
