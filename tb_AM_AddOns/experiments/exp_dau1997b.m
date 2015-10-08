function [masker BWtest] = exp_dau1997b(fs,nFig);
% function [masker BWtest] = exp_dau1997b(fs,nFig);
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 07/10/2015
% Last update on: 07/10/2015 
% Last use on   : 07/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

if nargin < 2
    nFig = 2; 
end

switch nFig
    case 2
    % To recreate Fig.3

        fc = 5000;
        dur = 10;
        SPL = 65;
        BWtest = [10 100 250 500 800 1000 2500 5000 10000];
        for i = 1:length(BWtest)
            masker(:,i) = AM_random_noise(fc-BWtest(i)/2,fc+BWtest(i)/2,SPL+3,dur,fs); 
        end
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
