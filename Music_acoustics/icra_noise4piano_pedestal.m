function outsig = icra_noise4piano_pedestal(insig,fs,RMSrel,SNR,fcmin,fcmax)
% function outsig = icra_noise4piano_pedestal(insig,fs,RMSrel,SNR,fcmin,fcmax)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 16/03/2016
% Last update on: 16/03/2016 
% Last use on   : 16/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    SNR = 0;
end

if nargin < 5
    fc = erbspacebw(80, 8000);
    fcmin = min(fc);
end

if nargin < 6
    fc = erbspacebw(80, 8000);
    fcmax = max(fc); 
end

rms_Buf = rmsdb_sec(insig,fs,350e-3);
rms_Buf = max(rms_Buf); % max peak value

erbmin = freqtoaud(min(fcmin))-0.5;
fmin   = audtofreq(erbmin);
erbmax = freqtoaud(max(fcmax))+0.5;
fmax   = audtofreq(erbmax);

noise  = AM_random_noise(fmin,fmax,60,5,fs); % 5-seconds white noise

pede_sample = auditoryfilterbank(noise,fs);
pede_sample = pede_sample(1:length(insig),:);
for i = 1:size(pede_sample,2)
    pede_sample(:,i) = From_dB(RMSrel(i))*pede_sample(:,i);
end
pede_sample = sum(pede_sample,2);
outsig = setdbspl(pede_sample,rms_Buf-SNR+100); % should be respect to the maximum value...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
