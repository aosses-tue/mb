function STI_and_SNRenv_DEMO
% function STI_and_SNRenv_DEMO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

bPlotSNRenv = 0; % if set to 1, the processing is very slow

load Danish_HINT_10sentence_samples_22kHz

nSents = 10; % number of concatenated sentences 
x_unprocessed = [];
for q = 1:nSents
    x_unprocessed = [ x_unprocessed sentenceArray{q}'];
    
end

fs = 22050;

x_unprocessed = resample(x_unprocessed,fs,22050);

% The input SNRs
inputSNRs = -6:2:6;

for n = 1:2
    % Load the noise file
    
    switch n
        case 1
            noise_name = 'SSN_HINT_22kHz';
        case 2
            noise_name = 'SSN_HINT_modReduced_-20dB.wav';
    end
    
    [noise_glob fs_n bits] = wavread(noise_name);

    if (fs_n~= fs)
        
        noise_glob = resample(noise_glob,fs ,fs_n);
        fs_n = fs;
    end
    
    
    if (fs ~= fs_n)
        error('sampling frequency not right')
    end
    x_org = x_unprocessed;
    
    Ts = 1/fs;
    T = length(x_org)/fs;
    t = 0:Ts:T;
    t = t(1:end-1);
    N = length(t);
    
    % scale the speech file to 65 dB SPL

    SPL = 65;
   
    speech= x_org/std(x_org)*10^((SPL)/20);
    
    noise_glob = noise_glob(1:N)';

    for k = 1:length(inputSNRs )
        
        noise = noise_glob;
        
        noise = noise/std(noise) * 10^((SPL-inputSNRs(k))/20)'; % scale the noise level
        
        if size(noise) ~= size(speech)
            noise = noise';
        end
        
        mixture = noise + speech;
        
        if bPlotSNRenv
            plotModExcitationPtn_for_Band =1;% plots the modulation excitation pattern in the 2-kHz band
            tmp(k,n)  = SNRenv_OctaveBand(mixture,noise,fs,plotModExcitationPtn_for_Band);
            SNRenvs(k,n) = tmp(k,n).SNRenv;
        end
        
        plotMTF_for_Band = 0;% plots the MTF in the 2-kHz band
        tmpSTIs(k,n) = sSTI(speech,mixture,fs,plotMTF_for_Band);
        STIs(k,n)= tmpSTIs(k,n).STI;
        
    end
     
end

figure
if bPlotSNRenv
    subplot(1,2,1)
end
plot(inputSNRs,STIs,'s', 'markersize',10)
xlabel('Input SNR (dB)')
ylabel('STI')
ylim([0 1])
xlim([-7 7])

if bPlotSNRenv
    subplot(1,2,2)
    plot(inputSNRs,SNRenvs,'o', 'markersize',10)
    ylabel('SNRenv (dB)')
    xlabel('Input SNR (dB)')
    xlim([-7 7])
    ylim([0 25])
end

legend('SSN','Modulation filtered SSN')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
