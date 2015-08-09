function y = r20150807_update_LNN_AMTControl
% function y = r20150807_update_LNN_AMTControl
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/08/2015
% Last update on: 06/08/2015 
% Last use on   : 07/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bEnergySum = 1;
bDoCreateLNN = 1;
bFFTusingPsySound = 1;

if bEnergySum
    
    N = 4096;
    K = N/2;
    fs = 44100; % K = 11025
    f = (1:K)/K * (fs/2); 
    SPLperHz = 40;
    BW = max(f);
    XdB = SPLperHz*ones(K,1);
    lvl = SPLperHz + 10*log10(BW);
    
    X = From_dB(XdB);
    lvl2b = rmsdb_freqdomain(f,X);
    
    df = f(2) - f(1);
    
    X2 = X.^2;
    lvl2 = 10*log10( df*sum(X2) );
    
    disp('')
    
end

if bDoCreateLNN
    Finf = 950; 
    Fsup = 1050; 
    BW = Fsup - Finf;
    fc = (Finf + Fsup)/2;
    Mdept = 0; % if 0 then Fmod is not important
    Fmod = 4;
    dur = 0.5; % [s]
    SPL = 70;
    fs = 44100;
    N_iterations = 10;
    % LNN_noise_BW(N_iterations,fc,BW,SPL,dur,fs,Fmod,Mdept); % to save waveforms
    [LNN y] = LNN_noise_BW(N_iterations,fc,BW,SPL+3,dur,fs,Fmod,Mdept);
    % sound(LNN,fs); 
    h = gcf;
    
    N = length(y);
    K = N/2;
    
    % windowtype = 'rectangular';
    windowtype = 'hanning';
    dBFS = 100;
    figure;
    freqfft2(y,K,fs,windowtype,dBFS);
     
    [z zdB f] = freqfft2(y,K,fs,windowtype,dBFS);
    lvl1 = rmsdb_freqdomain(f,z);
    
    [lnn lnndB] = freqfft2(LNN,K,fs,windowtype,dBFS);
    lvl2 = rmsdb_freqdomain(f,lnn);
    
    fprintf('Gaussian-noise: y = %.3f or %.3f\n',rmsdb(y)+dBFS,lvl1);
    fprintf('Low-noise noise: y = %.3f or %.3f\n',rmsdb(LNN)+dBFS,lvl2);
    
    figure; 
    plot(f,zdB,f,lnndB); grid on
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [dB SPL]')
    ylim([0 SPL+5])
    
    h(end+1) = gcf;
    
end

if bFFTusingPsySound
    f1 = 'D:\Output\LNNnoise-fc-1000_BW-100_Fmod-4_Mdept-0_SPL-70.wav';
    f2 = 'D:\Output\randomnoise-fc-1000_BW-100_Fmod-4_Mdept-0_SPL-70.wav';
    PS_CL_01_20150312_reduced(f1,f2);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
