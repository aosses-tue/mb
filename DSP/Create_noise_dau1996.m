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
% Last update on: 23/10/2014 % Update this date manually
% Last use on   : 23/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dB_SPL_noise    = options.dB_SPL_noise;
fs              = options.fs;
[t_silence_bef, t_duration, t_silence_aft] = Create_noise_dau1996_default(nTag);
t_tmp           = 200e-3;
Ntmp            = round(t_tmp*fs);

Nnoise          = round(fs*(t_duration+t_tmp));

y   = wgn(Nnoise,1,1);

y   =  y(:); % ensures it is a column vector

Wn = [20 5000]/(fs/2); % Normalised cutoff frequency        
[b,a] = butter(2,Wn); % 8th-order
tmp_y = filtfilt(b,a,y); % Linear-phase implementation

y = y(Ntmp+1:end);
tmp_y = tmp_y(Ntmp+1:end);

% TO DEBUG:
% tt = (1:length(y))/fs; figure; plot(tt,y,tt,tmp_y); xlim([-0.1 0.2]), legend('bef filt','after filt')
%
% TO DEBUG, applying a cosine ramp:
% yw = y; Ramp = rampup(5e-3*fs); L = length(Ramp); yw(1:L) = yw(1:L).*Ramp; tmp_yw = filtfilt(b,a,yw);
% figure; plot(tt,y,tt,tmp_yw); xlim([-0.1 0.2]), legend('bef filt','after filt with ramp')

y   = setdbspl(tmp_y,dB_SPL_noise);
innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end

options = Ensure_field(options,'bSave_noise',0);

if nargout == 0 | options.bSave_noise 
    Wavwrite(innoise,fs,filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
