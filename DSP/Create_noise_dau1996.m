function [innoise, fs] = Create_noise_dau1996(nTag,filename,options)
% function [innoise, fs] = Create_noise_dau1996(nTag,filename,options)
%
% 1. Description:
%
% 2. Stand-alone example:
%       nTag = 3; % 600-ms white noise, no onset; no offset
%       filename = 'test';
%       options.dB_SPL_noise = 77;
%       options.fs = 48000;
%       innoise = Create_noise_dau1996(nTag,filename,options);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 16/03/2015 % Update this date manually
% Last use on   : 16/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SPL_noise       = options.dB_SPL_noise;
fs              = options.fs;
[t_silence_bef, t_duration, t_silence_aft] = Create_noise_dau1996_default(nTag);
t_tmp           = 200e-3;
Ntmp            = round(t_tmp*fs);

dur             = t_silence_bef + t_duration + t_silence_aft;
Nnoise          = round(fs*(t_duration+t_tmp));

Finf    = 20;
Fsup    = 5000;

method  = 0; 

switch method
    case 0 % default: since 16/03/2015
        
        Fmod    = 4; % not important if Mdept is 0
        Mdept   = 0;
        dBFS    = 100;
        SPL     = SPL_noise+3;
        attack  = 0; % no cos_ramp
        release = 0; % no cos_ramp
        y = AM_random_noise(Finf,Fsup,SPL,dur,fs,Fmod,Mdept,dBFS);
        y = cos_ramp(length(y), fs, attack, release);
        
    case 2  % this implementation introduced an emphasis either at the onset
            % or at the offset of the white noise... (check it out)
        y   = wgn(Nnoise,1,1);
        y   =  y(:); % ensures it is a column vector
        Wn = [Finf Fsup]/(fs/2); % Normalised cutoff frequency
        [b,a] = butter(2,Wn); % 8th-order
        tmp_y = filtfilt(b,a,y); % Linear-phase implementation
        y = y(Ntmp+1:end);
        tmp_y = tmp_y(Ntmp+1:end);
        y   = setdbspl(tmp_y,SPL_noise);
        attack = 25;
        release = 25;
        y = cos_ramp(length(y), fs, attack, release);
end

innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end

options = Ensure_field(options,'bSave_noise',0);

if nargout == 0 | options.bSave_noise 
    Wavwrite(innoise,fs,filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
