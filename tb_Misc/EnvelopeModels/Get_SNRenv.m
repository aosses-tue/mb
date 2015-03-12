function h = Get_SNRenv(insig,noise_name,fs,testSNRs)
% function h = Get_SNRenv(insig,noise_name,fs,testSNRs)
%
% File analysed step by step by Alejandro Osses
% Adapted from: STI_and_SNRenv_DEMO.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    testSNRs = -6:2:6; % The input SNRs
end

insig = resample(insig,fs,22050);

for n = 1
    % Load the noise file
    [noise fs_n bits] = wavread(noise_name);

    if (fs_n~= fs)
        noise = resample(noise,fs ,fs_n);
        fs_n = fs;
    end
    
    if (fs ~= fs_n)
        error('sampling frequency not right')
    end
    
    Ts = 1/fs;
    T = length(insig)/fs;
    t = 0:Ts:T;
    t = t(1:end-1);
    N = length(t);
    
    % scale the speech file to 65 dB SPL

    SPL = 65;
   
    speech  = insig/std(insig)*10^((SPL)/20);
    noise   = noise(1:N)';

    for k = 1:length(testSNRs )
        
        noise = noise/std(noise) * 10^((SPL-testSNRs(k))/20)'; % scale the noise level
        
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
plot(testSNRs,STIs,'s', 'markersize',10)
xlabel('Input SNR (dB)')
ylabel('STI')
ylim([0 1])
xlim([-7 7])

if bPlotSNRenv
    subplot(1,2,2)
    plot(testSNRs,SNRenvs,'o', 'markersize',10)
    ylabel('SNRenv (dB)')
    xlabel('Input SNR (dB)')
    xlim([-7 7])
    ylim([0 25])
end

h = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
