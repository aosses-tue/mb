function [STIs MTF_mean MTF_std outs] = Get_STI(insig,noise_name,fs,testSNRs,SPL)
% function [STIs MTF_mean MTF_std outs] = Get_STI(insig,noise_name,fs,testSNRs,SPL)
%
% File analysed step by step by Alejandro Osses
% Adapted from: STI_and_SNRenv_DEMO.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    testSNRs = -6:2:6; % The input SNRs
end

if nargin < 5
    SPL     = 65; % The input SNRs
end

fsanalysis  = 22050;
insig       = resample(insig,fsanalysis,fs);
fs       	= fsanalysis;

[noise fs_n bits]   = Wavread(noise_name); % Load the noise file

if fs_n~= fs
    
    noise   = resample(noise,fs,fs_n);
    fs_n    = fs;
    
end

Ts  = 1/fs;
T   = length(insig)/fs;
t   = 0:Ts:T;
t   = t(1:end-1);
N   = length(t);

% scale the speech file to SPL in dB
speech  = insig/std(insig)*10^((SPL)/20);
noise   = noise(1:N)';

MTF_mean   = [];
MTF_std    = [];

for k = 1:length( testSNRs )

    noise = noise/std(noise) * 10^((SPL-testSNRs(k))/20)'; % scale the noise level

    if size(noise) ~= size(speech)
        noise = noise';
    end

    mixture = noise + speech;

    plotMTF_for_Band = 0;% plots the MTF in the 2-kHz band
    tmpSTIs(k,1)= sSTI(speech,mixture,fs,plotMTF_for_Band);
    STIs(k,1)   = tmpSTIs(k,1).STI;
    MTF_mean(k,1:14)= mean( transpose( tmpSTIs(k,1).MTF ) );
    MTF_std(k,1:14) = std( transpose( tmpSTIs(k,1).MTF ) );

end
  
outs.modBands = tmpSTIs(k,1).modBands; 
outs.audBands = tmpSTIs(k,1).audBands; 

if nargout == 0
    figure
    plot(testSNRs,STIs,'s', 'markersize',10)
    xlabel('Input SNR (dB)')
    ylabel('STI')
    ylim([0 1])
    xlim([-7 7])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
