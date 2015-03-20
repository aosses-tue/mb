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
bPart2 = 0;
bPart3 = 1;

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

stim_phase = [10 20 40 80 160 320];

figure;
errorbar(stim_phase, mean([85+SRTs1; 85+SRTs2]),std([85+SRTs1; 85+SRTs2]), 'o-');
xlabel('Signal duration [ms]')
ylabel('Masked threshold [dB SPL]')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\';
dirResults = [dir 'Results' delim];
dirStimuli = [dir 'Stimuli' delim];

% file3 = [dir 'dau1996b-I2-3AFC-results-AO-part1.apr.xml'];
file1   = [dirResults 'dau1996b-I2-3AFC-multi-1-of-2-results.apr.xml']; 
file2   = [dirResults 'dau1996b-I2-3AFC-multi-2-of-2-results.apr.xml'];
file3   = [dirResults 'dau1996b-I2-3AFC-multi-1-of-2-results-2.apr.xml']; 
file4   = [dirResults 'dau1996b-I2-3AFC-multi-2-of-2-results-2.apr.xml'];
file5   = [dirResults 'dau1996b-I2-3AFC-multi-1-of-2-results-3.apr.xml']; 
file6   = [dirResults 'dau1996b-I2-3AFC-multi-2-of-2-results-3.apr.xml'];

stim1     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-0-pi.wav'];
stim2     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-5_10-pi.wav'];
stim3     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-10_10-pi.wav'];
stim4     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-13_10-pi.wav'];
stim5     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-15_10-pi.wav'];
stim6     = [dirStimuli 'dau1996b_expI2_stim-5ms-85-phase-20_10-pi.wav'];

stimnoise = [dir 'Stimuli' delim 'dau1996b_expI_noisemasker.wav'];

x1          = Wavread(stim1); x1(end) = []; % it has one more sample
x2          = Wavread(stim2); x2(end) = []; 
x3          = Wavread(stim3); x3(end) = []; 
x4          = Wavread(stim4); x4(end) = []; 
x5          = Wavread(stim5); x5(end) = []; 
x6          = Wavread(stim6); x6(end) = []; 

[xnoise fs] = Wavread(stimnoise);
t = ( 0:length(xnoise)-1 )/fs;

x_minmax = minmax(xnoise');
tonset = 115e-3;
sonset = round(tonset*fs);

Norm_factor = 1/max(abs(x_minmax));
PhaseN = Norm_factor* xnoise(sonset) ;

offset_y = 0.25;
M = max(offset_y*(6+1),1);

figure;
plot(t, Norm_factor*xnoise - 0.4); hold on, grid on
plot(t, 0.1*Norm_factor*x1 + offset_y*6,'r')
plot(t, 0.1*Norm_factor*x2 + offset_y*5,'r')
plot(t, 0.1*Norm_factor*x3 + offset_y*4,'r')
plot(t, 0.1*Norm_factor*x4 + offset_y*3,'r')
plot(t, 0.1*Norm_factor*x5 + offset_y*2,'r')
plot(t, 0.1*Norm_factor*x6 + offset_y*1,'r')
legend('Frozen noise','1-kHz tone diff phase');

plot([tonset tonset]            ,[-1 M],'k'   ,'LineWidth',2)
plot([tonset+5e-3 tonset+5e-3]  ,[-1 M],'k--' ,'LineWidth',2)
xlim([tonset-2.5e-3 tonset+7.5e-3])
ylim([-1 M])

title(sprintf('Stimuli. Noise with phase of %.2f x pi at signal onset.',PhaseN))
xlabel('Time [s]')
ylabel('Normalised amplitude')

opts.N4mean = 10;
opts.bPlot  = 0;
dBoffset    = 85;
% [SRTs1 SRTdev1] = quick_staircases(file1,opts);
% [SRTs2 SRTdev2] = quick_staircases(file2,opts);
out1 = quick_staircases_multi(file1,opts); SRTs1 = out1.SRT+dBoffset; SRTdev1 = out1.Stdev;
out2 = quick_staircases_multi(file2,opts); SRTs2 = out2.SRT+dBoffset; SRTdev2 = out2.Stdev;
out3 = quick_staircases_multi(file3,opts); SRTs3 = out3.SRT+dBoffset; SRTdev3 = out3.Stdev;
out4 = quick_staircases_multi(file4,opts); SRTs4 = out4.SRT+dBoffset; SRTdev4 = out4.Stdev;
out5 = quick_staircases_multi(file5,opts); SRTs5 = out5.SRT+dBoffset; SRTdev5 = out5.Stdev;
out6 = quick_staircases_multi(file6,opts); SRTs6 = out6.SRT+dBoffset; SRTdev6 = out6.Stdev;

stim_phase = [0 1/2 1 1.3 3/2 2];

MeanV = mean([SRTs1 SRTs2; ...
              SRTs3 SRTs4; ...
              SRTs5 SRTs6])
StdV = std([SRTdev1 SRTdev2; ...
            SRTdev3 SRTdev4; ...
            SRTdev5 SRTdev6 ])
figure;
errorbar(stim_phase, MeanV,StdV, 'o-');
hold on
plot(stim_phase,[SRTs1 SRTs2],'rx','LineWidth',2)
plot(stim_phase,[SRTs3 SRTs4],'rx','LineWidth',2)
plot(stim_phase,[SRTs5 SRTs6],'rx','LineWidth',2)

xlabel('Signal phase x \pi [rad]')
ylabel('Masked threshold [dB SPL]')
grid on
xlim([0 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
