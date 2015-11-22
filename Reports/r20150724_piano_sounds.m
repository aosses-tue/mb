function r20150724_piano_sounds
% function r20150724_piano_sounds
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/07/2015
% Last update on: 19/11/2015 
% Last use on   : 19/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bDoWaveforms = 0;
bDoSTFT_f0  = 1;
bDoSTFT     = 0;
bDoMFB = 1;

bPrepareFigures = 0;

dir = [Get_TUe_paths('Databases') 'dir01-Instruments' delim 'Piano' delim '04-PAPA' delim];
f1 = [dir 'F1' delim 'NS19-F1.wav']; % Used for Expose
f2 = [dir 'F1' delim 'JBS36-F1.wav']; % Used for Expose
% dir = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\F1-44100-Hz\';
% f1 = [dir 'NS19-F1-noisered.wav']; % Only first 5 seconds
% f2 = [dir 'JBS36-F1-noisered.wav']; % Only first 5 seconds

f1 = [dir 'A4' delim 'JBS73-A4.wav']; 
f2 = [dir 'A4' delim 'NS19-A4.wav']; 

tmax = 3; % s
fmax = 2000; % Hz

sens = 50e-3;
G    = 5;
Cal  = 1/(G*sens);  % 1 =  94 dB
Cal  = Cal/2;       % 1 = 100 dB 

title1 = 'Type NS19';
title2 = 'Type JBS36';

[x1 fs1] = Wavread(f1);
[x2 fs2] = Wavread(f2);

x1 = Cal  *x1;
x2 = Cal/2*x2;

fs = fs1;
t1 = ( 1:length(x1) )/fs;
t2 = ( 1:length(x2) )/fs;

if bDoWaveforms
    figure;
    %subplot(1,2,1)
    plot(t1,2*x1), grid on
    xlim([0 tmax])
    title( sprintf('%s',title1) )
    xlabel('Time [s]') 
    ylabel('Pressure [Pa]')
    h(1) = gcf;
    
    figure;
    plot(t2,2*x2), grid on
    xlim([0 tmax])
    title( sprintf('%s',title2) )
    xlabel('Time [s]') 
    ylabel('Pressure [Pa]')
    h(2) = gcf;
end

if bDoSTFT
    
    nfft = 4096*2;
    wlen = nfft/2;
    overlap = 75;
    nwtype = 4; % Hamming window
    
    figure
    stft(x1, fs, nfft, wlen, overlap, nwtype);
    xlabel('Time [s]') 
    title( sprintf('%s',title1) )
    ylim([0 fmax])
    xlim([0 tmax])
    h(3) = gcf;
    
    figure;
    stft(x2, fs, nfft, wlen, overlap, nwtype);
    xlabel('Time [s]') 
    title( sprintf('%s',title2) )
    ylim([0 fmax])
    xlim([0 tmax])
    h(4) = gcf;
    
end

if bPrepareFigures
    hM.I_TitleInAxis = 0;
    hM.I_KeepColor = 0;
    hM.I_FontSize = 18;
    hM.I_Width = 8;
    hM.I_Height = 4;
    
    if bDoSTFT
        Figure2paperfigureT(h(1:2),1,2,hM) % manually stored
    end
    if bDoSTFT
        Figure2paperfigureT(h(3:4),1,2,hM) % manually stored
    end
end

if bDoMFB
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options.fc2plot = 90;
    % options.fc2plot = 800;
    % options.fc2plot = 1500;
    options.Title = title1;
    insig1 = x1(1:tmax*fs);
    demo_dau1997b_one_sound(insig1,fs,options);
    
    options.Title = title2;
    insig2 = x2(1:tmax*fs);
    demo_dau1997b_one_sound(insig2,fs,options);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
