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
% Last update on: 11/02/2015 % Update this date manually
% Last use on   : 11/02/2015 % Update this date manually
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

if bAnalyseModelFiles
    
    dir = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\';
    
    % Ac. mode 2, far
    files = {   [dir 'mode-1-v_1.wav'],
                [dir 'mode-1-v_1-ane.wav'],
                [dir 'mode-1-v_1-rev.wav']};
	[xtot fs] = Wavread(files{1});
    xane = Wavread(files{2});
    xrev = Wavread(files{3});
    xsum = xane + xrev;
    
    t = ( 1:length(xtot) )/fs;
    
    ha = [];
    figure;
    subplot(4,1,1)
    plot(t,xtot);
    ha(end+1) = gca;
    title('Total, Ac.mode 2, distant mic.')
    
    subplot(4,1,2)
    plot(t,xane);
    ha(end+1) = gca;
    title('Direct')
    
    subplot(4,1,3)
    plot(t,xrev);
    ha(end+1) = gca;
    title('Reflected')
    
    subplot(4,1,4)
    plot(t, abs(xsum-xtot) );
    ha(end+1) = gca;
    title('Error')
    
    linkaxes(ha,'x')
    
    xlim([0 4*0.607])
    
    % Ac. mode 2, close mic
    files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2.wav',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2-ane.wav',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-1-v_2-rev.wav'};
	[xtot fs] = Wavread(files{1});
    xane = Wavread(files{2});
    xrev = Wavread(files{3});
    xsum = xane + xrev;
    
    t = ( 1:length(xtot) )/fs;
    
    ha = [];
    figure;
    subplot(4,1,1)
    plot(t,xtot);
    ha(end+1) = gca;
    title('Total, Ac.mode 2, close mic')
    
    subplot(4,1,2)
    plot(t,xane);
    ha(end+1) = gca;
    title('Direct')
    
    subplot(4,1,3)
    plot(t,xrev);
    ha(end+1) = gca;
    title('Reflected')
    
    subplot(4,1,4)
    plot(t, abs(xsum-xtot) );
    ha(end+1) = gca;
    title('Error')
    
    linkaxes(ha,'x')
    
    xlim([0 4*0.607])
    
    % ac-mode 4
    files = {   'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1.wav',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1-ane.wav',
                'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test\Stage2\mode-3-v_1-rev.wav'};
    
    [xtot fs] = Wavread(files{1});
    xane = Wavread(files{2});
    xrev = Wavread(files{3});
    xsum = xane + xrev;
    
    t = ( 1:length(xtot) )/fs;
    
    ha = [];
    figure;
    subplot(4,1,1)
    plot(t,xtot);
    ha(end+1) = gca;
    title('Total, Ac.mode 4')
    
    subplot(4,1,2)
    plot(t,xane);
    ha(end+1) = gca;
    title('Direct')
    
    subplot(4,1,3)
    plot(t,xrev);
    ha(end+1) = gca;
    title('Reflected')
    
    subplot(4,1,4)
    plot(t, abs(xsum-xtot) );
    ha(end+1) = gca;
    title('Error')
    
    linkaxes(ha,'x')
    
    xlim([0 4*0.3])
    
    % Ac. mode 4, close
    files = {   [dir 'mode-3-v_2.wav'],
                [dir 'mode-3-v_2-ane.wav'],
                [dir 'mode-3-v_2-rev.wav']};
	[xtot fs] = Wavread(files{1});
    xane = Wavread(files{2});
    xrev = Wavread(files{3});
    xsum = xane + xrev;
    
    t = ( 1:length(xtot) )/fs;
    
    ha = [];
    figure;
    subplot(4,1,1)
    plot(t,xtot);
    ha(end+1) = gca;
    title('Total, Ac.mode 4, close mic.')
    
    subplot(4,1,2)
    plot(t,xane);
    ha(end+1) = gca;
    title('Direct')
    
    subplot(4,1,3)
    plot(t,xrev);
    ha(end+1) = gca;
    title('Reflected')
    
    subplot(4,1,4)
    plot(t, abs(xsum-xtot) );
    ha(end+1) = gca;
    title('Error')
    
    linkaxes(ha,'x')
    
    xlim([0 4*0.3])
    
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
