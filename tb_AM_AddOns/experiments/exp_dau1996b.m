function [masker,insig] = exp_dau1996b(fs)
% function [masker insig] = exp_dau1996b(fs)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 23/09/2015
% Last update on: 23/09/2015 
% Last use on   : 23/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

f = 1000;
durramp = 5; % ms
lvl = 75; % arbitrary

nExperiment = 3; % signal integration

if nExperiment == 3
    
    masker = AM_random_noise(20,5000,60,600e-3,fs); 
    masker = Do_cos_ramp(masker,fs,durramp);
    masker = setdbspl(masker,77);
    
    dur = [10 20 40]; % ms
    
    insig = [];
    for i = 1:3
        insigtmp = Create_sin(f,dur(i)*1e-3,fs);
        insigtmp = Do_cos_ramp(insigtmp,fs,durramp);
        insigtmp = setdbspl(insigtmp,lvl);
        insigtmp = [Gen_silence(100e-3,fs); insigtmp; Gen_silence((500-dur(i))*1e-3,fs)]; % 600 ms 
        insig = [insig insigtmp];
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
