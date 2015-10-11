function [masker,insig] = exp_jepsen2008(fs,nFig)
% function [masker insig] = exp_jepsen2008(fs,nFig)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fs      = 44100;
%       nFig    = 3.2;
%       [masker insig] = exp_jepsen2008(fs,nFig);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/10/2015
% Last update on: 05/10/2015 
% Last use on   : 05/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

f = 1000;
durramp = 5; % ms
lvl = 75; % arbitrary

if nargin < 2
    nFig = 3.2; % signal integration
end

switch nFig
    case 3.2 
    % To recreate Fig.3, right panel: intensity discrimination task

        dur = 10; % buffer of 10 s
        maskerbuf = AM_random_noise(100,10000,60+3,dur,fs); 
        
        dur = 0.5;
        rampupdn = 50; % ms
        masker = AM_random_noise(100,10000,60+3,dur,fs); 
        masker = Do_cos_ramp(masker,fs,rampupdn);
        
        if nargout == 0
            f1 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-60-dB-dur-500-ms.wav'];
            f2 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-42-dB-dur-500-ms.wav'];
            f3 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-60-dB-dur-10-s.wav'];
            f4 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-42-dB-dur-10-s.wav'];
            Wavwrite(       masker,fs,f1); % sound(masker,fs);
            Wavwrite(gaindb(masker,-18),fs,f2); % sound(gaindb(masker,-18),fs);
            Wavwrite(       maskerbuf,fs,f3); % sound(masker,fs);
            Wavwrite(gaindb(maskerbuf,-18),fs,f4); % sound(gaindb(masker,-18),fs);
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
