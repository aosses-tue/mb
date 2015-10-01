function r20151002_update_dau
% function r20151002_update_dau
%
% 1. Description:
%       This script is a copy of r20150320_update_dauM.m, but re-run in Sept
%       2015, where some parts of the code were optimised (and updated).
% 
% 2. Stand-alone example:
%       r20151002_update_dau;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150320_update, r20150320_updateN
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/09/2015
% Last update on: 28/09/2015 
% Last use on   : 29/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

bDiary = 0;
Diary(mfilename,bDiary);

bPart1          = 1;
bPart2          = 1; % Simulation
bPart2recalculate = 1;
if bPart1
    bCreateAudio= 0; % run this first
end

opts.calc_method = 2;

bPlot = 1;

h               = [];

tmp = Get_TUe_paths('outputs');
pathaudiosrc  = [tmp 'new' delim '48000' delim];
pathaudiodst1 = [tmp 'new' delim '44100' delim];% 'D:\Output\r20141031_update_dau_et_al20150318\';
pathaudiodst2 = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\'; % APEX results
pathaudiodst3 = [tmp 'new' delim 'Experiment' delim];
output_dir = [Get_TUe_paths('lx_Text') 'lx2015-03-20-update' delim 'Figures' delim];

test_onset_ref  = [-20:5:10]*1e-3;

x = [];
offset = 0.2;

if bPart1
    % figure;
    for i = 1:10
        fname   = [pathaudiosrc 'dau1996b_expIB0_stim-10ms-76-' num2str(i) '.wav'];
        fname2  = sprintf('%sdau1996b_expIB0_stim-10ms-76-onset-%.0f-ms',pathaudiodst1,test_onset_ref(i)*1000+50);
        [tmp fs]= Wavread(fname);
        x(:,i) = tmp;
        t = (0:length(tmp)-1)/fs;

        % plot(t*1000,tmp+offset*(i-1)), hold on

        y = resample(tmp,44100,fs);
        if bCreateAudio
            Wavwrite(y,44100,fname2);
        end
    end

    N = size(y,1);
    fname = [pathaudiodst1 'dau1996b_expI_noisemasker.wav']; % existing already
    [tmp fs] = Wavread(fname);
    M = size(tmp,1);
    x(:,i+1) = [tmp; tmp(1:N-M-1)];
else
    % fill in filenames
    % generates 50ms of silence prior to the noise and truncates to noise length
end

dirResults = [pathaudiodst2 'Results' delim];
dirStimuli = [pathaudiodst3 delim]; % same as [pathaudiodst2 'Stimuli' delim];

file1   = [dirResults 'dau1996b-IIC-3AFC-multi-1-of-2-results-3.apr.xml']; 
file2   = [dirResults 'dau1996b-IIC-3AFC-multi-1-of-2-results-4.apr.xml'];
file3   = [dirResults 'dau1996b-IIC-3AFC-multi-2-of-2-results-1.apr.xml']; 
file4   = [dirResults 'dau1996b-IIC-3AFC-multi-2-of-2-results-2.apr.xml'];

stim1   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-30-ms.wav']; % -20
stim2   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-35-ms.wav']; % -15
stim3   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-40-ms.wav']; % -10
stim4   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-45-ms.wav']; % -5
stim5   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav']; % 0
stim6   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-55-ms.wav']; % 5
stim7   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-60-ms.wav']; % 10
stimnoise = [dirStimuli 'dau1996b_expI_noisemasker.wav'];

filenames{1} = stim1;
filenames{2} = stim2;
filenames{3} = stim3;
filenames{4} = stim4;
filenames{5} = stim5;
filenames{6} = stim6;
filenames{7} = stim7;
filenames{8} = stimnoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2
    
    if bPart2recalculate
        opts.bExpIC0 = 1; % Backward masking

        opts.method = 'dau1996';
        opts.pathaudio = pathaudiodst3;
        opts.filenames = filenames;
        r20141031_update_dau_et_alN(opts);

        opts.method = 'dau1996a';
        opts.pathaudio = pathaudiodst3;
        opts.filenames = filenames;
        r20141031_update_dau_et_alN(opts);

    else

        error('Not updated yet')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Backward masking:
        % path = 'D:\Output\r20141031_update_dau_et_al20150318\';
        % f1 = [path 'mue-2015-03-18-at-19h-59m-40s.mat']; % no limitation
        % f2 = [path 'mue-2015-03-18-at-20h-11m-24s.mat']; % limitation
        % 
        % load(f1); % loads mue
        % 
        % dB_SPL          = [10 16:10:76]; % Get this from processing
        % test_onset_ref  = [-100:40:-20,-15:5:15]*1e-3; % Get this from processing
        % 
        % criterion_corr = 8.5;% options.criterion_corr;
        % 
        % Thres   = demo_dau1996b_decision(mue(3:end,:),dB_SPL,criterion_corr);
        % 
        % load(f2); % loads mue
        % ThresO  = demo_dau1996b_decision(mue(3:end,:),dB_SPL,criterion_corr);
        % 
        % figure;
        % plot(test_onset_ref(3:end)*1000,Thres,'rx-'), grid on, hold on
        % plot(test_onset_ref(3:end)*1000,ThresO,'bo-','MarkerFaceColor','b');
        % 
        % errorbar(stim_onset,MeanV,StdV,'k-')
        % legend('no limit','limit = 10','Measured Thresholds','Location','NorthWest')
        % xlabel('Signal onset relative to masker onset [ms]')
        % ylabel('Masked thresholds [dB]')
        % 
        % h(end+1) = gcf;
    end
end

Save_all_figures(h,output_dir);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
