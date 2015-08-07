function y = r20150807_update_LNN_AMTControl
% function y = r20150807_update_LNN_AMTControl
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/08/2015
% Last update on: 06/08/2015 
% Last use on   : 06/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bDoCreateLNN = 1;

if bDoCreateLNN
    Finf = 50; %950;
    Fsup = 5000; % 1050;
    BW = Fsup - Finf;
    fc = (Finf + Fsup)/2;
    Mdept = 0; % if 0 then Fmod is not important
    Fmod = 4;
    dur = 0.5; % [s]
    SPL = 70;
    fs = 44100;
    N_iterations = 10;
    % LNN_noise_BW(N_iterations,fc,BW,SPL,dur,fs,Fmod,Mdept); % to save waveforms
    [LNN y] = LNN_noise_BW(N_iterations,fc,BW,SPL,dur,fs,Fmod,Mdept);
    % sound(LNN,fs); 
    
    K = length(y)/2;
    h = Get_window('hamming',2*K);
    corr = length(h)*mean(h);
    corrWindow_dB = log10(corr) + 0.272;
%     h = Get_window('rectangular',2*K);
%     corr = 1;
    opts.fs = fs;
    [z zdB f] = freqfft(h.*y,K,opts);
    % [zi zdBi] = freqfft(yinterim,K,opts);
    [lnn lnndB] = freqfft(h.*LNN,K,opts);
    
    RMSf_corr = rmsdb(z)-10*log10(length(y))+corrWindow_dB;
    RMSf_corrLNN = rmsdb(lnn)-10*log10(length(LNN))+corrWindow_dB;
    
    fprintf('Gaussian-noise: y = %.3f or %.3f\n',rmsdb(y),RMSf_corr);
    fprintf('Low-noise noise: y = %.3f or %.3f\n',rmsdb(LNN),RMSf_corrLNN);
    
    figure; 
    plot(f,zdB + corrWindow_dB,f,lnndB + corrWindow_dB); grid on
    
    disp('')
    
%     
%     factor2 = 10^(rmsdb(y2)/10);
%     z2m = factor2 * z2;
%     ZRMS = rmsdb(z2m) - 10*log10(length(y2));
    
end



% 
% W=1.438, V=-18.5
% W=1.438, V=-19.2
% W=1.439, V=-19.7
% W=1.439, V=-20.1
% y = -33.000 or -37.125
% r20150807_update_LNN_AMTControl
% Diary (log) for function/script r20150807_update_LNN_AMTControl disabled. Change bDoDiary to 1 to enable this feature...
% W=2.819, V=-5.9 (for no iteration)
% W=1.652, V=-5.9
% W=1.521, V=-10.8
% W=1.475, V=-13.0
% W=1.456, V=-14.8
% W=1.447, V=-16.2
% W=1.444, V=-17.3
% W=1.442, V=-18.1
% W=1.441, V=-18.7
% W=1.440, V=-19.1
% W=1.441, V=-19.5
% y = -33.000 or -37.435
% r20150807_update_LNN_AMTControl
% Diary (log) for function/script r20150807_update_LNN_AMTControl disabled. Change bDoDiary to 1 to enable this feature...
% W=2.799, V=-6.0 (for no iteration)
% W=1.646, V=-6.0
% W=1.517, V=-11.0
% W=1.477, V=-13.3
% W=1.457, V=-15.1
% W=1.447, V=-16.6
% W=1.443, V=-17.8
% W=1.443, V=-18.7
% W=1.443, V=-19.3
% W=1.443, V=-19.7
% W=1.443, V=-20.1
% y = -33.000 or -37.031
% r20150807_update_LNN_AMTControl
% Diary (log) for function/script r20150807_update_LNN_AMTControl disabled. Change bDoDiary to 1 to enable this feature...
% W=2.765, V=-6.0 (for no iteration)
% W=1.629, V=-6.0
% W=1.506, V=-11.0
% W=1.470, V=-13.3
% W=1.456, V=-15.2
% W=1.449, V=-16.6
% W=1.446, V=-17.7
% W=1.444, V=-18.6
% W=1.443, V=-19.1
% W=1.442, V=-19.5
% W=1.441, V=-19.8
% y = -33.000 or -36.985
% edit Roughness_offline.m
% r20150807_update_LNN_AMTControl
% Diary (log) for function/script r20150807_update_LNN_AMTControl disabled. Change bDoDiary to 1 to enable this feature...
% W=2.746, V=-6.0 (for no iteration)
% W=1.641, V=-6.0
% W=1.521, V=-11.0
% W=1.480, V=-13.2
% W=1.461, V=-15.0
% W=1.451, V=-16.5
% W=1.447, V=-17.6
% W=1.445, V=-18.4
% W=1.444, V=-19.0
% W=1.444, V=-19.4
% W=1.444, V=-19.8
% y = -33.000 or -37.315
% log10(corrdB)
% 
% ans =
% 
%     4.0424
% 
% corr = length(h)*mean(h);
%     
%     RMSf_corr = rmsdb(z)-10*log10(length(y))-log10(corr);
% fprintf('y = %.3f or %.3f\n',rmsdb(y),RMSf_corr)
% y = -33.000 or -41.357
% RMSf_corr = rmsdb(z)-10*log10(length(y))+log10(corr);
% fprintf('y = %.3f or %.3f\n',rmsdb(y),RMSf_corr);
% y = -33.000 or -33.272
% RMSf_corr = rmsdb(z)-10*log10(length(y))+log10(corr)+0.272;
% r20150807_update_LNN_AMTControl
% Diary (log) for function/script r20150807_update_LNN_AMTControl disabled. Change bDoDiary to 1 to enable this feature...
% W=2.851, V=-5.9 (for no iteration)
% W=1.643, V=-5.9
% W=1.517, V=-10.9
% W=1.476, V=-13.1
% W=1.456, V=-14.9
% W=1.446, V=-16.3
% W=1.442, V=-17.4
% W=1.440, V=-18.2
% W=1.440, V=-18.9
% W=1.440, V=-19.3
% W=1.440, V=-19.7
% y = -33.000 or 3.396
% corr_additional = -10*log10(length(y))+corr;
%     RMSf_corr = rmsdb(z) + corr_additional;
% fprintf('y = %.3f or %.3f\n',rmsdb(y),RMSf_corr);
% y = -33.000 or 3.396
% corr = 10*log10( length(h)*mean(h) )+0.272;
% RMSf_corr = rmsdb(z)-10*log10(length(y))+corr;
% corr = 10*log10( length(h)*mean(h);
%     
%     RMSf_corr = rmsdb(z)-10*log10(length(y))+log10(corr)+0.272;



f1 = 'D:\Output\LNNnoise-fc-1000_BW-100_Fmod-4_Mdept-0_SPL-70.wav';
f2 = 'D:\Output\randomnoise-fc-1000_BW-100_Fmod-4_Mdept-0_SPL-70.wav';
PS_CL_01_20150312_reduced(f1,f2)

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
