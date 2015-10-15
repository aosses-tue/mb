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
%       fs      = 44100;
%       nFig    = 7;
%       [masker insig] = exp_jepsen2008(fs,nFig);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/10/2015
% Last update on: 12/10/2015 
% Last use on   : 12/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

f = 1000;
durramp = 5; % ms

if nargin < 2
    % nFig = 3.2; % signal integration
    nFig = 7; % forward-masking
end

switch nFig
    case 3.2 
    % To recreate Fig.3, right panel: intensity discrimination task

        dur = 10; % buffer of 10 s
        maskerbuf = AM_random_noise(100,10000,60+3,dur,fs); 
        
        for i = 1 % 1:4
            dur = 0.5;
            rampupdn = 50; % ms
            masker = AM_random_noise(100,10000,60+3,dur,fs); 
            masker = Do_cos_ramp(masker,fs,rampupdn);

            if nargout == 0
                f1 = sprintf('%sjepsen2008-BW-at-60-dB-dur-500-ms-%.0f.wav',Get_TUe_paths('outputs'),i);
                f2 = sprintf('%sjepsen2008-BW-at-42-dB-dur-500-ms-%.0f.wav',Get_TUe_paths('outputs'),i);
                f3 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-60-dB-dur-10-s.wav'];
                f4 = [Get_TUe_paths('outputs') 'jepsen2008-BW-at-42-dB-dur-10-s.wav'];
                Wavwrite(       masker,fs,f1); % sound(masker,fs);
                Wavwrite(gaindb(masker,-18),fs,f2); % sound(gaindb(masker,-18),fs);
                Wavwrite(       maskerbuf,fs,f3); % sound(masker,fs);
                Wavwrite(gaindb(maskerbuf,-18),fs,f4); % sound(gaindb(masker,-18),fs);
            end
        end
        
    case 7
        dur = 10; % for buffer
        fmasker1 = 4000;
        fmasker2 = 2400;
        ftarget = 4000;
        durtarget = 10e-3;
        sil1 = 500e-3;
        sil2 = sil1 + 30e-3;
        durtargettot = 1;
        sil1c = durtargettot - sil1 - durtarget;
        sil2c = durtargettot - sil2 - durtarget;
        
        lvl = 60; % arbitrary
        wintarget = 1; % 1 = hanning
        ymasker1 = Create_sin(fmasker1,dur,fs);
        ymasker2 = Create_sin(fmasker2,dur,fs);
        ymasker1 = setdbspl(ymasker1,lvl);
        ymasker2 = setdbspl(ymasker2,lvl);
        
        ytarget = Create_sin(ftarget,durtarget,fs,wintarget);
        ytarget = setdbspl(ytarget,lvl);
        ytarget1 = [Gen_silence(sil1,fs); ytarget; Gen_silence(sil1c,fs)];
        ytarget2 = [Gen_silence(sil2,fs); ytarget; Gen_silence(sil2c,fs)];
        
        f1 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-60-dB-dur-250-ms-onset-%.0f-ms.wav',Get_TUe_paths('outputs'),ftarget,sil1*1000);
        f2 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-60-dB-dur-250-ms-onset-%.0f-ms.wav',Get_TUe_paths('outputs'),ftarget,sil2*1000);
        f3 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-masker-60-dB-dur-10-s.wav',Get_TUe_paths('outputs'),fmasker1);
        f4 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-masker-60-dB-dur-10-s.wav',Get_TUe_paths('outputs'),fmasker2);
        
        if nargout == 0
            Wavwrite(ytarget1,fs,f1); 
            Wavwrite(ytarget2,fs,f2); 
            Wavwrite(ymasker1,fs,f3); 
            Wavwrite(ymasker2,fs,f4); 
        end
    
        % on-frequency:
        lvls = 30:80;
        
        % off-frequency
        lvls = [60 70 80 85];
        
    case 99
        dur = 10; % for buffer
        fmasker1 = 800;
        ftarget  = 800;
        durtarget= 10e-3;
        sil1     = 200e-3;
        durtargettot = 0.7;
        sil1c    = durtargettot - sil1 - durtarget;
        
        lvl      = 60; % arbitrary
        wintarget=  1; % 1 = hanning
        ymasker1 = Create_sin(fmasker1,dur,fs);
        ymasker1 = setdbspl(ymasker1,lvl);
        
        ytarget  = Create_sin(ftarget,durtarget,fs,wintarget);
        ytarget  = setdbspl(ytarget,lvl);
        ytarget1 = [Gen_silence(sil1,fs); ytarget; Gen_silence(sil1c,fs)];
                
        f1 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-60-dB-dur-250-ms-onset-%.0f-ms.wav',Get_TUe_paths('outputs'),ftarget,sil1*1000);
        f2 = sprintf('%sjepsen2008-fig7-%.0f-Hz-tone-masker-60-dB-dur-10-s.wav',Get_TUe_paths('outputs'),fmasker1);
        
        if nargout == 0
            Wavwrite(ytarget1,fs,f1); 
            Wavwrite(ymasker1,fs,f2); 
        end
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
