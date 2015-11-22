function r20151119_piano_sounds
% function r20151119_piano_sounds
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

bDoSTFT_f0  = 0;

bDoWaveforms = 0;
bDoSTFT     = 0;
bPrepareFigures = 0;

notetmp.note = 'F'; notetmp.octave = 3;
% notetmp.note = 'A'; notetmp.octave = 4;
f0target = note2freq(notetmp);
note2test = [notetmp.note num2str(notetmp.octave)];
dir = [Get_TUe_paths('Databases') 'dir01-Instruments' delim 'Piano' delim '04-PAPA' delim note2test delim];

opts.bExtension = 0;
files = Get_filenames(dir,'wav',opts);

tmax = 3; % s
fmax = 2000; % Hz

sens = 50e-3;
G    = 5;
Cal  = 1/(G*sens);  % 1 =  94 dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoSTFT_f0
    
    for i = 1:length(files)
        title1 = files{i};
        fullfile{i} = [dir files{i}];

        [insig fs] = Wavread([fullfile{i} '.wav']);
    
        nfft = 4096*4;
        wlen = nfft/2;
        overlap = 75;
        nwtype = 4; % Hamming window
        [y1 f t1] = stft(insig, fs, nfft, wlen, overlap, nwtype);

        idxt1 = 1:length(t1); % find(t1>0.3 & t1<1.4);
        idxf  = find(f>0 & f<1.3*f0target); %find(f>0.7*f0target & f<1.3*f0target);

        y1dB = To_dB(abs(y1));
        [xx idx] = max(y1dB(idxf,:));
        fmax1 = f(idx(idxt1));

        figure; 
        plot(t1(idxt1),fmax1); grid on
        title(title1)
        
        disp('')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(files)
    title1 = files{i};
    fullfile{i} = [dir files{i}];
    [insig fs] = Wavread([fullfile{i} '.wav']);
    
    f0(i) = il_get_f0([fullfile{i} '.txt']);

end

dir_new = sprintf('%snorm-%.0f-Hz%s',dir,f0target,delim);
Mkdir(dir_new);

for i = 1:length(files)
    
    [insig fs] = Wavread([fullfile{i} '.wav']);
    
    Perc = round(100*(f0target-f0(i))/f0target); 
    outsig = Do_pitch_stretch(insig,fs,Perc,'percentage');
    
    fullfile_new{i} = [dir_new files{i} '.wav'];
    Wavwrite(outsig,fs,fullfile_new{i});
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% if bPrepareFigures
%     hM.I_TitleInAxis = 0;
%     hM.I_KeepColor = 0;
%     hM.I_FontSize = 18;
%     hM.I_Width = 8;
%     hM.I_Height = 4;
%     
%     if bDoSTFT
%         Figure2paperfigureT(h(1:2),1,2,hM) % manually stored
%     end
%     if bDoSTFT
%         Figure2paperfigureT(h(3:4),1,2,hM) % manually stored
%     end
% end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function f0 = il_get_f0(file)

fileID = fopen(file);
C = textscan(fileID,'%.1f %s');
fclose(fileID);
f0 = C{1,1};
