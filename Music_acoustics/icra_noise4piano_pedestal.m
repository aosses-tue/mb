function [outsig signal_plus_pede] = icra_noise4piano_pedestal(insig,fs,RMSrel,SNR,fcmin,fcmax)
% function [outsig signal_plus_pede] = icra_noise4piano_pedestal(insig,fs,RMSrel,SNR,fcmin,fcmax)
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

method = 1; % method 0 is respect to max RMS
            % method 1 loudness balancing respect to ICRA noise (SNR = 0 dB) then gain is applied to match the target SNR

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

rms_Buf = rmsdb_sec(insig,fs,10e-3);
rms_Buf = max(rms_Buf); % max peak value

erbmin = freqtoaud(min(fcmin))-0.5;
fmin   = audtofreq(erbmin);
erbmax = freqtoaud(max(fcmax))+0.5;
fmax   = audtofreq(erbmax);

noise  = AM_random_noise(fmin,fmax,60,5,fs); % 5-seconds white noise

pede_sample = auditoryfilterbank(noise,fs);
% pede_sample = pede_sample(1:length(insig),:);
for i = 1:size(pede_sample,2)
    pede_sample(:,i) = From_dB(RMSrel(i))*pede_sample(:,i);
end
pede_sample = sum(pede_sample,2);

switch method
    case 0
        outsig = setdbspl(pede_sample,rms_Buf-SNR+100); % should be respect to the maximum value...
    case 1
        pede_sample = setdbspl(pede_sample,rms_Buf+100);
        [pede_sample gain lou] = il_loudness_balance(insig,pede_sample,fs); % 'output at 0 dB SNR'
        outsig = From_dB(-SNR) * pede_sample;
end

% pede_sample = pede_sample(1:length(insig),:);

signal_plus_pede = insig + outsig(1:length(insig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outsig gain lou] = il_loudness_balance(insig_ref,insig,fs)
% function [outsig gain lou] = il_loudness_balance(insig_ref,insig,fs)

dur_ana = 350e-3;
N = round(dur_ana * fs);
[inst_l, short_l] = Get_Loudness_MGB(insig_ref(1:N),fs);
loudness_max = max(short_l);

disp('Loudness balancing of pedestal noise (respect to ICRA sounds)')
if loudness_max < 1 | loudness_max > 25
    fprintf('\n\tAssumes that the piano onset occurs during the first %.1f [s] of the input signal\n\n',dur_ana);
end
    
gain = 0; % current gain dB
step = 1.5; % 1.5 dB
    
[inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
lou = max(short_loud);
if lou > loudness_max
    go_down = 1;
else
    go_down = 0;
end
       
switch go_down
    case 1
        while (lou > loudness_max & abs(lou-loudness_max) > 0.5)

            gain = gain - step;
            [inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
            lou = max(short_loud);

        end
    case 0
        while ( lou < loudness_max & abs(lou-loudness_max) > 0.5 )

            gain = gain + step;
            [inst_loud,short_loud] = Get_Loudness_MGB(From_dB(gain)*insig(1:round(dur_ana*fs)), fs);
            lou = max(short_loud);

        end

end
outsig = From_dB(gain)*insig;


