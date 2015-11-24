function y = r20151127_piano_sounds_noise
% function y = r20151127_piano_sounds_noise
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/11/2015
% Last update on: 24/11/2015 
% Last use on   : 24/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

bDoResampling   = 0;
bDoPianoSounds  = 1;
bDoSTFT     = 0;
bDoMFB = 1;

bPrepareFigures = 0;

dir = [Get_TUe_paths('Databases') 'dir01-Instruments' delim 'Piano' delim '04-PAPA' delim];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoResampling
    
    notetmp = { 'Dsh', 1; ...
                'F'  , 1; ... 
                'C'  , 2; ...
                'Ash', 2; ...
                'F'  , 3;...
                'C'  , 4;... 
                'A'  , 4;...
                'Csh', 5;... 
                'C'  , 6;... 
                'G'  , 6}; 
	fstarget = 44100;
    
    for i = 1:length(notetmp)
        
        ntmp.note   = notetmp{i,1};
        ntmp.octave = notetmp{i,2};
        
        noteS    = [ntmp.note num2str(ntmp.octave)];
        f0target = note2freq(ntmp);
                
        dirtarget{i} = sprintf('%s%s%snorm-%.0f-Hz%s',dir,noteS,delim,f0target,delim); % resampled-at-44100-Hz
        Resample2fs(dirtarget{i},fstarget);
        
        dstfolder = [dir '01-Tuned-at-44100-Hz-new' delim noteS delim];
        movefile([dirtarget{i} 'resampled-at-44100-Hz' delim],dstfolder);
        
        disp('')
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir44100 = [dir '01-Tuned-at-44100-Hz' delim];
f1 = [dir44100 'A4' delim 'GH05-A4.wav']; 

sens = 50e-3;
G    = 5;
Cal  = 1; %1/(G*sens);  % 1 =  94 dB
% Cal  = Cal/2;       % 1 = 100 dB (AMT convention)

if bDoPianoSounds
    
    [x1 fs] = Wavread(f1);
    x1 = Cal*x1;
    
    Create_piano_noise(x1,fs); % get plot
    [env noise w] = Create_piano_noise(x1,fs);
    
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

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
