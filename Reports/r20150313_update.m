function r20150313_update
% function r20150313_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2015
% Created on    : 12/03/2015
% Last update on: 12/03/2015 % Update this date manually
% Last use on   : 12/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bPart1 = 0;
bPart2 = 1;

close all

if bPart1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path    = Get_TUe_paths('db_voice_of_dragon');
dirm    = [path '02-Wav-files'           delim '2015-02-wav-files' delim '02-calibrated' delim];
dirp    = [path '03-Wav-files-predicted' delim '2015-02-wav-files' delim '02-calibrated' delim];
fnoise1 = [Get_TUe_paths('db_calfiles') 'track_03.wav'];

% f1 = [dirm 'meas-ac-4-dist-rev.wav']; label1 = 'meas-rev';
% % f1 = [dirp 'model-ac-4-dist-ane.wav']; label1 = 'mod-ane';
% f2 = [dirp 'model-ac-4-dist-rev.wav']; label2 = 'mod-rev';

f1      = '~/Documenten/Databases/dir03-Speech/dutch/list/alle-zinnen/wdz48.wav';
fnoise1 = '~/Documenten/Databases/dir03-Speech/dutch/list/alle-zinnen/wivineruis.wav';
fnoise2 = '~/Documenten/Databases/dir03-Speech/dutch/list/alle-zinnen/whitenoise-LISTf.wav';
label1  = 'SSN';
label2  = 'White';

PS_CL_10_20150312(f1,fnoise1);
% PS_CL_10_20150312(f2,fnoise1);
PS_CL_10_20150312(f1,fnoise2);

[insig1 fs]  = Wavread(f1);
[insig2 fs]  = Wavread(f1);
testSNRs     = -5:5:20;

SPL1 = rmsdb(insig1)+100;
SPL2 = SPL1; % rmsdb(insig2)+100;

[STI1 MTF_mean1 MTF_std1 out1] = Get_STI(insig1,fnoise1,fs,testSNRs,SPL1);
[STI2 MTF_mean2 MTF_std2]      = Get_STI(insig1,fnoise2,fs,testSNRs,SPL2);
modBands = out1.modBands;

figure
plot(testSNRs,[STI1 STI2],'s', 'markersize',10)
xlabel('Input SNR (dB)')
ylabel('STI')
ylim([0 1])
% xlim([-7 7])
legend(label1, label2);
grid on
% STI increases as the noise decreases (i.e, when SNR increase)

%% Answer 5: Plot MTF for the different input SNRs

figure
for i = 1:7
    subplot(3,3,i)
    errorbar(testSNRs-0.1,[MTF_mean1(:,i)],[MTF_std1(:,i)]), hold on
    errorbar(testSNRs+0.1,[MTF_mean2(:,i)],[MTF_std2(:,i)],'r')
    title(sprintf('Mod-band %.2f [Hz]',modBands(i)))
    grid on
end
legend(label1, label2)%,'Location','EastOutside');

figure;
for i = 8:14
    subplot(3,3,i-7)
    errorbar(testSNRs-0.1,[MTF_mean1(:,i)],[MTF_std1(:,i)]), hold on
    errorbar(testSNRs+0.1,[MTF_mean2(:,i)],[MTF_std2(:,i)],'r')
    title(sprintf('Mod-band %.2f [Hz]',modBands(i)))
    grid on
end
legend(label1, label2)%,'Location','EastOutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\';
file1 = [dir 'dau1996b-3AFC-multi-results-AO.apr.xml'];
file2 = [dir 'dau1996b-3AFC-multi-results-RK.apr.xml'];

stim1 = [dir 'Stimuli' delim 'dau1996b_expI3_stim01-85.wav']; % Hanning (check this out)

[x1 fs] = Wavread(stim1);
t1 = ( 0:length(x1)-1 )/fs;

figure;
plot(t1,x1);

opts.N4mean = 10;
opts.bPlot = 0;
[SRTs1 SRTdev1] = quick_staircases(file1,opts);
[SRTs2 SRTdev2] = quick_staircases(file2,opts);

stim_durations = [10 20 40 80 160 320];

figure;
errorbar(stim_durations, mean([85+SRTs1; 85+SRTs2]),std([85+SRTs1; 85+SRTs2]), 'o-');
xlabel('Signal duration [ms]')
ylabel('Masked threshold [dB SPL]')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
