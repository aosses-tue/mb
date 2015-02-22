function r20150213_update
% function r20150213_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 11/02/2015
% Last update on: 13/02/2015 % Update this date manually
% Last use on   : 13/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

%% These parameters were taken from r20150130_update:

bCountPeriods = 0;
bCreateModelFiles = 0;
bAnalyseModelFiles = 1;

dir_audio = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\1-referentie-BP-filtered\';
% dir_audio_mod_ane   = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data\at44100Hz\';
% dir_audio_modelled  = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-reflections\at44100Hz\';
% output_dir  = [Get_TUe_paths('outputs') 'Hummer-audio-files' delim];
    
ac_mode     = [2 4]; % [2 3 4 5];
ac_mode_idx = 1:length(ac_mode);
cal_level   = [60 78]; % [60 75 78 82]; % ac mode 2, 3, 5
delta_dB_far = [-6 0]; % distance distant mic = 2 * close mic
ver_anechoic = [2 3]; %[2 1 3 1];     % ac-mode-5, take 1 = contaminated with ac mode 4
ver_reverberant = [1 1]; % [1 2 1 3]; % ac-mode-5

field = 3;

%%
if bCountPeriods == 1
    for i = field

        % anechoic
        for j = ac_mode_idx

            files = sprintf('modus %.0f_v%.0f-%.0f.wav',ac_mode(j)-1,ver_anechoic(ac_mode_idx(j)), i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio files];

            [o1 o2] = Count_times_above_thr(fullfile);
            rotA(j) = o1;
        end

        % Reverberant conditions:
        for j = ac_mode_idx

            files = sprintf('modus %.0f_v%.0f-%.0f.wav',ac_mode(j)-1,ver_reverberant(ac_mode_idx(j)), i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio files];

            [o1 o2] = Count_times_above_thr(fullfile);
            rotR(j) = o1; 

        end
    end
end

if bCreateModelFiles
    
    % ac-mode 2
    files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_1.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_1-ane.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_1-rev.txt'};
	for i = 1:length(files)
        VoD_one_predicted_file_from_txt(files{i});
    end
    
    % ac-mode 2, close
    files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2-ane.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2-rev.txt'};
	for i = 1:length(files)
        VoD_one_predicted_file_from_txt(files{i});
    end
    
    % ac-mode 4
    files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1-ane.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1-rev.txt'};
       
	for i = 1:length(files)
        VoD_one_predicted_file_from_txt(files{i});
    end
     
	% ac-mode 4, close
	files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_2.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_2-ane.txt',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_2-rev.txt'};
	for i = 1:length(files)
        VoD_one_predicted_file_from_txt(files{i});
    end
            
end

T = [0 0.6 0 0.3 0];
if bAnalyseModelFiles
    
    % dir = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\';
    
    % Added p, 13/02/2015:
    dir = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test2\';
    
    h = [];
    for i = [2 4]
        % Ac. mode 2, far
        files = {   [dir 'Stage4\mode-' num2str(i-1) '-v_1-rev.wav'], ...
                    [dir 'Stage3\mode-' num2str(i-1) '-v_1-s1r.wav'], ...
                    [dir 'Stage3\mode-' num2str(i-1) '-v_1-s1m.wav'], ...
                    [dir 'Stage3\mode-' num2str(i-1) '-v_1-s2r.wav'], ...
                    [dir 'Stage3\mode-' num2str(i-1) '-v_1-s2m.wav']};

        %[dir 'Stage4\mode-1-v_1.wav'],

        [xtot fs] = Wavread(files{1});

        x1r = Wavread(files{2});
        x1m = Wavread(files{3});
        x2r = Wavread(files{4});
        x2m = Wavread(files{5});

        xane = x1r + x2r;
        xrev = x1m + x2m;

        xsum = xane + xrev;

        rmstot = rmsdb(xsum);
        rmsane = rmsdb(xane);
        rmsrev = rmsdb(xrev);

        t = ( 1:length(xtot) )/fs;

        ha = [];
        figure;
        subplot(3,1,1)
        plot(t,xtot);
        ha(end+1) = gca;
        title(sprintf('Ac.mode %.0f, distant mic., reverberant condition',i))

        subplot(3,1,2)
        plot(t,xane);
        ha(end+1) = gca;
        title(sprintf('Direct sound, rel.level = %.2f dB',rmsane-rmstot));
        ylabel('Amplitude')
        
        subplot(3,1,3)
        plot(t,xrev);
        ha(end+1) = gca;
        title('Reflected')
        title(sprintf('Reflected sound, rel.level = %.2f dB',rmsrev-rmstot));

        xlabel('Time (s)')
        % subplot(4,1,4)
        % plot(t, abs(xsum-xtot) );
        % ha(end+1) = gca;
        % title('Error')

        linkaxes(ha,'x')

        xlim([0 4*T(i)])

        h(end+1) = Figure2paperfigure(gcf);

    end
    
    acmode = [2 4];
    for i = 1:length(h)
        Saveas(h(i),sprintf('waveforms-ac-mode-%.0f',acmode(i)));
    end
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
