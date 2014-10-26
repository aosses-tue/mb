function [innoise, fs] = Create_noise_dau1996(nTag,filename,options)
% function [innoise, fs] = Create_noise_dau1996(nTag,filename,options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       nTag = 3; % 600-ms white noise, no onset; no offset
%       filename = 'test';
%       options.dB_SPL_noise = 77;
%       options.fs = 48000;
%       innoise = Create_noise_dau1996(nTag,filename,options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 23/10/2014 % Update this date manually
% Last use on   : 23/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dB_SPL_noise    = options.dB_SPL_noise;
fs              = options.fs;
[t_silence_bef, t_duration, t_silence_aft] = Create_noise_dau1996_default(nTag);

Nnoise      = round(fs*t_duration);

y = wgn(Nnoise,1,1);

y   =  y(:); % ensures it is a column vector

Wn = [20 5000]/(fs/2); % Normalised cutoff frequency        
[b,a] = butter(4,Wn); % 8th-order
y = filtfilt(b,a,y); % Linear-phase implementation

y   = setdbspl(y,dB_SPL_noise);
innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end

options = Ensure_field(options,'bSave_noise',0);

if nargout == 0 | options.bSave_noise 
    Wavwrite(innoise,fs,filename);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
