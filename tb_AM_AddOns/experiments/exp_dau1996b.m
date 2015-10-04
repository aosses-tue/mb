function [masker,insig] = exp_dau1996b(fs,nFig)
% function [masker insig] = exp_dau1996b(fs,nFig)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fs      = 44100;
%       nFig    = 14;
%       [masker insig] = exp_dau1996b(fs,nFig);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 23/09/2015
% Last update on: 23/09/2015 
% Last use on   : 02/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    fs = 44100;
end

f = 1000;
durramp = 5; % ms
lvl = 75; % arbitrary

if nargin < 2
    nFig = 3; % signal integration
end

switch nFig
    case {3,4}
    % To recreate Fig.3 and Fig. 4: signal integration experiment

        masker = AM_random_noise(20,5000,60,600e-3,fs); 
        masker = Do_cos_ramp(masker,fs,durramp);
        masker = setdbspl(masker,77);

        dur = [10 20 40 70 150]; % ms

        insig = [];
        for i = 1:length(dur)
            insigtmp = Create_sin(f,dur(i)*1e-3,fs);
            insigtmp = Do_cos_ramp(insigtmp,fs,durramp);
            insigtmp = setdbspl(insigtmp,lvl);
            insigtmp = [Gen_silence(100e-3,fs); insigtmp; Gen_silence((500-dur(i))*1e-3,fs)]; % 600 ms 
            insig = [insig insigtmp];
        end
        
    case {14,15}
        
        masker = AM_random_noise(20,5000,60,5,fs); % 5-seconds buffer
        masker = setdbspl(masker,77);

        dur = 5e-3; % s
        tonset = [0 10 20 50 100]*1e-3;
        durn = 300*1e-3; %
        durs = durn - dur;
        durnN = durn*fs;
        
        insig = [];
        for i = 1:length(tonset)
            insigtmp = Create_sin(f,dur,fs,'hanning'); % hanning-windowed over entire duration
            insigtmp = setdbspl(insigtmp,lvl);
            insigtmp = [Gen_silence(tonset(i),fs); insigtmp; Gen_silence((durn-tonset(i)-dur),fs)]; 
            if length(insigtmp) ~= durnN
                insigtmp = [insigtmp; 0];
            end
            insig = [insig insigtmp];
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
