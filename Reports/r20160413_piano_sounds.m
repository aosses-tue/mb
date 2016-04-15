function r20160413_piano_sounds
% function r20160413_piano_sounds
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
% Original file name: r20151119_piano_sounds.m
% Created on    : 13/04/2016
% Last update on: 13/04/2016 
% Last use on   : 13/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bSave = 1;

bDo_f0_only_check = 0;
bDo_f0_shift = 1;
method_f0 = 1;

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

opts.bExtension = 0;

%%%
sens = 50e-3;   % Given by Antoine
G    = 5;       % Given by Antoine

Cal  = 1/(G*sens);  % 1 =  94 dB
Cal  = Cal/2;       % 1 = 100 dB (AMT convention)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDo_f0_only_check
    
    % dir     = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Tuned-at-44100-Hz' delim note2test delim];
    % dir     = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Tuned-at-44100-Hz' delim note2test delim];
    dir     = [Get_TUe_data_paths('piano') '04-PAPA' delim '01-Sounds' delim note2test delim 'norm-' num2str(round(f0target)) '-Hz-ESPRIT' delim];
    files   = Get_filenames(dir,['wav']);
    for i = 1:length(files)
        title1 = files{i};
        fullfile{i} = [dir files{i}];
        [insig fs] = Wavread([fullfile{i}]);

        tolerance = 10;
        N   = round(0.1*fs); 
        [outsig Fi Ai] = il_get_ESPRIT_piano(insig,fs,N,f0target,tolerance);

    end
end

if bDo_f0_shift
    dir     = [Get_TUe_data_paths('piano') '04-PAPA' delim '01-Sounds' delim note2test delim];
    files   = Get_filenames(dir,'wav',opts);

    for i = 1:length(files)
        title1 = files{i};
        fullfile{i} = [dir files{i}];
        try
            [insig fs] = Wavread([fullfile{i} '.wav']);
        catch
            [insig fs] = Wavread([fullfile{i}]);
        end

        switch method_f0
            case 0
                f0(i) = il_get_f0([fullfile{i} '.txt']);
            case 1

            N   = round(0.1*fs); % N has to be greater than 5*L, arbitrarily chosen
            [outsig Fi Ai] = il_get_ESPRIT_piano(insig,fs,N,f0target);
            f0(i) = Fi(1);
            
            fprintf('\tf0 = %.3f [Hz] (Sound %s)\n\n',f0(i),files{i});
        end
    end

    switch method_f0 
        case 0
            dir_new = sprintf('%snorm-%.0f-Hz%s',dir,f0target,delim);
        case 1
            dir_new = sprintf('%snorm-%.0f-Hz-ESPRIT%s',dir,f0target,delim);
    end
    Mkdir(dir_new);

    for i = 1:length(files)

        try
            [insig fs] = Wavread([fullfile{i} '.wav']);
        catch
            [insig fs] = Wavread([fullfile{i}]);
        end
        
        switch method_f0
            case 0
                Perc = round(100*(f0target-f0(i))/f0(i)); % shift to be defined respect to the initial f0 
            case 1
                Perc = 100*(f0target-f0(i))/f0(i); 
        end
            
        %[outsig f0out]  = Do_pitch_stretch(insig,Perc,'percentage',fs);
        
        if Perc > 3
            outsig = insig;
            while Perc > 3
                Perc = 3;
                outsig  = Do_pitch_stretch(outsig,Perc,'percentage',fs);
                tolerance = 5;
                [xx Fi Ai f0_new] = il_get_ESPRIT_piano(outsig,fs,N,f0target,tolerance);
                Perc = 100*(f0target-f0_new)/f0_new; 
                % Perc = Perc_still;
            end
            insig = outsig;
        end
        
        if Perc > 0.5
            [outsig f0out]  = Do_pitch_stretch(insig,Perc,'percentage',fs);
        else
            outsig = insig;
            disp([files{i} ' untouched'])
        end
        
        fullfile_new{i} = [dir_new files{i} '.wav'];
        if bSave
            Wavwrite(outsig,fs,fullfile_new{i});
        end

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

function [outsig Fi Ai f0] = il_get_ESPRIT_piano(insig,fs,N,f0target,tolerance)

if nargin < 5
    tolerance = 10;
end

[insig_max Ni] = max(abs(insig)); % max of the waveform (manually computed)
Ni = Ni + round(0.020*fs); % 20-ms after the max
Nf  = Ni+N-1; 
insig_t = insig(Ni:Nf);
p = 180;
[outsig Fi Ai] = Get_ESPRIT_analysis(insig_t,p,N,fs);

[Amax idxmax] = max(Ai); % first 10 partials
idxFi = find( 100*abs(Fi/f0target-1)<tolerance ); % 10 %
factor2div = 1;

partial_nr = 2;

while length(idxFi) == 0
    idxFi = find( 100*abs((Fi/partial_nr)/f0target-1)<tolerance/partial_nr ); % second harmonic
    factor2div = partial_nr;
    partial_nr = partial_nr + 1;
end

Ait = Ai(idxFi);
Fit = Fi(idxFi)/factor2div; % if factor2divide ~= 0, then 'virtual pitch' (when using ESPRIT) 
[Amax idxmax] = max(Ait); % first 10 partials
f0 = Fit(idxmax);

if factor2div == 1
    fprintf('\tf0 = %.3f [Hz]\n\n',f0);
else
    fprintf('\tf0 = %.3f [Hz] - ''virtual pitch'' \n\n',f0);
end