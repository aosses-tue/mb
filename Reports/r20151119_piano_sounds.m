function r20151119_piano_sounds
% function r20151119_piano_sounds
%
% 1. Description:
%       bDoSTFT_f0 - checks out f0 for every piano sound of the refister set
%                    by notetmp
%       bDo_f0_shift - does the shift of the piano sounds. It is required 
%                    that the txt files with the f0 notes (as recorded has
%                    been already generated (by 20/11/2015 done for A4 and F3
%                    registers).
% 
% 2. Stand-alone example:
%       r20151119_piano_sounds;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20150724_piano_sounds.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 19/11/2015
% Last update on: 22/11/2015 
% Last use on   : 22/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bDoSTFT_f0  = 0;
bDo_f0_shift = 0;
bDoEnv = 1;

% notetmp.note = 'Dsh'; notetmp.octave = 1; % Not yet % very difficult to label
% notetmp.note = 'F'; notetmp.octave = 1; 
% notetmp.note = 'C'; notetmp.octave = 2;
% notetmp.note = 'Ash'; notetmp.octave = 2;
% notetmp.note = 'F'; notetmp.octave = 3;
notetmp.note = 'C'; notetmp.octave = 4; 
% notetmp.note = 'A'; notetmp.octave = 4;
% notetmp.note = 'Csh'; notetmp.octave = 5; 
% notetmp.note = 'C'; notetmp.octave = 6; 
% notetmp.note = 'G'; notetmp.octave = 6; 

f0target = note2freq(notetmp);
note2test = [notetmp.note num2str(notetmp.octave)];
dir = [Get_TUe_paths('Databases') 'dir01-Instruments' delim 'Piano' delim '04-PAPA' delim note2test delim];

opts.bExtension = 0;
files = Get_filenames(dir,'wav',opts);
% files = Get_filenames([dir 'norm-117-Hz' delim],'wav',opts);

%%%
sens = 50e-3;   % Given by Antoine
G    = 5;       % Given by Antoine

Cal  = 1/(G*sens);  % 1 =  94 dB
Cal  = Cal/2;       % 1 = 100 dB (AMT convention)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoSTFT_f0
    
    % This process was not automated. Run it for every piano sound and write
    % down in a txt file with the same name of the wav file the fundamental
    % frequency read from the plots (e.g. GH05-A4.txt if GH05-A4.wav was plotted
    
    for i = 1:length(files)
        title1 = files{i};
        fullfile{i} = [dir files{i}];
        
        [insig fs] = Wavread([fullfile{i} '.wav']);  % non-calibrated audio files
        
        nfft = 4096*4;
        wlen = nfft/2;
        overlap = 75;
        nwtype = 4; % Hamming window
        [y1 f t1] = stft(insig, fs, nfft, wlen, overlap, nwtype);

        idxt1 = 1:length(t1); % find(t1>0.3 & t1<1.4);
        idxf  = find(f>0.8*f0target & f<1.3*f0target); %find(f>0.7*f0target & f<1.3*f0target);

        y1dB = To_dB(abs(y1));
        [xx idx] = max(y1dB(idxf,:));
        idx = idx + min(idxf)-1;
        fmax1 = f(idx(idxt1));

        figure; 
        plot(t1(idxt1),fmax1); grid on
        title(title1)
        
        disp('')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDo_f0_shift
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDoEnv
    for i = 1:length(files)
        title1 = files{i}; 
        % fullfile{i} = [dir files{1}]; % using first file (before adjustment)
        fullfile{i} = sprintf('%snorm-%.0f-Hz%s%s',dir,f0target,delim,files{1}); % using first file (after adjustment)
        [insig fs] = Wavread([fullfile{i} '.wav']);
        
        method = 3;
        switch method
            case 1 % Hilbert
                y = abs(hilbert( insig ));
                
            case 2
                
                y = abs(hilbert( insig ));
                % timeconstant = 1/(2*pi*f0);
                % f0 = 1/(2*pi*timeconstant);
                f0 = f0target*0.1;
                [b a] = IRIfolp(f0,fs);
                y = filter(b,a,y);
                
            case 3
                
                yin = abs(hilbert( insig ));
                [b, a] = butter(4,20/(fs/2),'low');
                
                % % A look at the group delay introduced by the filter shows that the delay is frequency-dependent.
                % grpdelay(b,a,4096,fs); % Large group delay around cut-off frequency of the butter filter

                ylp = filtfilt(b,a,yin);
                
                % Aslow = 8; Apeak = 6; % Parameters used by RK
                Aslow = 2; Apeak = 0; 
                
                ydiff = yin - Aslow*ylp;  % Makes Aslow*ylp well below 0 (except for the onset)
                                        % figure; plot(ydiff)
                ydiff = max(ydiff,0);
                
                yout = ydiff*Apeak;
                
                y = yin + yout;
                                
                % figure; freqz(b,a,4096);
                figure;
                subplot(2,1,1)
                plot(yin); hold on
                plot(yout,'r');
                ha = gca;
                % sound(ytot,fs)
                       
        end
        
        % idxstart = 900;
        % y2plot = ylp(idxstart:end)+yout(1:end-idxstart+1);
        y2plot = ylp+yout;
        subplot(2,1,2);
        plot(insig(1:end)); hold on
        plot(y2plot(1:end),'r');
        ha(end+1)=gca;
        linkaxes(ha,'x')
        xlim([2.2e5 2.8e5])
        disp('')
    end
     
end

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
