function meas_20131025_PR_stimuli_analysis
% function meas_20131025_PR_stimuli_analysis
%
% 1. Description:
%   It processes sequences of tones centred at 104, 110, 117, 124 and 131 Hz.
%   Four sequences are read, 2 for F0mod, 2 for ACE (wav files have same F0 
%   but different roving values) 
% 
%   freq is mapped into freq_warped. freq_warped is linearly interpolated
%   between fi and fs of each FAT-channel
% 
% Programmed by Alejandro Osses
% Created in: 2013
% Last updated on: 21/05/2014
% Last used on   : 21/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bSave = 1;
if isunix
    directory = ['/home/alejandro/Documenten/Meas/Meas/Electrodograms/20131026-checking-JL-stimuli/'];
else
    directory = ['D:\SVN-KU-Leuven\alejandro\Meas\Meas\Electrodograms\20131026-checking-JL-stimuli\'];
end

addpath(directory)
addpath([directory 'PR_Stimuli-JL/'])

load([directory 'tg_p-q-structures_ACE_Simulink-var.mat']);

seqACE      = read_aseq([directory 'Sequence-104Hz-JL-ACE-highest.aseq']); % Sequence with tones at 104, 110, 117, 124 and 131 Hz
seqF0m      = read_aseq([directory 'Sequence-104Hz-JL-F0m-highest.aseq']);
seqACE_low  = read_aseq([directory 'Sequence-104Hz-JL-ACE-lowest.aseq']);
seqF0m_low  = read_aseq([directory 'Sequence-104Hz-JL-F0m-lowest.aseq']);

seqACEpatient = MapElectrodogram2PatientBoudaries(seqACE,p);
seqF0mpatient = MapElectrodogram2PatientBoudaries(seqF0m,p);
seqACEpatient_low = MapElectrodogram2PatientBoudaries(seqACE_low,p);
seqF0mpatient_low = MapElectrodogram2PatientBoudaries(seqF0m_low,p);

H_dB        = [];
H_dB_F0m    = [];
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_ACE_104_Hz.wav');
[~          , H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_ACE_110_Hz.wav');
[~          , H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_ACE_117_Hz.wav');
[~          , H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_ACE_124_Hz.wav');
[~          , H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_ACE_131_Hz.wav');

[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_104_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_110_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_117_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_124_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_131_Hz.wav');

xaxis_tones     = [0 650 1380 2100 2830];
PlotSize        = [99 274 1811 425];
hFig    = [];
hAx     = [];
Plot_sequence(seqACEpatient)
title('ACE, balanced stimuli for JL + 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hACE_high = gca;
 
Plot_sequence(seqACEpatient_low)
title('ACE, balanced stimuli for JL - 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hACE_low = gca;

Plot_sequence(seqF0mpatient)
title('F0mod, balanced stimuli for JL + 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hF0m_high = gca;

Plot_sequence(seqF0mpatient_low)
title('F0mod, balanced stimuli for JL - 5 dB')
hFig(end+1)     = gcf;
hAx(end+1)      = gca;
hF0m_low        = gca;

for i = 1:length(hFig)
    set(hFig(i),'Position', PlotSize)
end


hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
dest_folder     = [hoofd_folder 'ci_pooled/Figures/'];

linkaxes([hACE_high hF0m_high hACE_low hF0m_low])
ylim([0 20])

for j = 1:length(xaxis_tones)
    
    for i = 1:length(hFig)
        if bSave
            if mod(i,2) == 1
                prefix = 'h';
            else
                prefix = 'l';
            end
            
            figure(hFig(i));
            xlim([xaxis_tones(j)-50 xaxis_tones(j)+600])
            rectangle('Position',[xaxis_tones(j)-50 0.01 50 20-0.05],'FaceColor','w'), hold on
            if i <=2 % ACE plots
                plot(xaxis_tones(j)-50+H_dB(:,j)*50,freq_warped,'r')
                legend('seq','FFT')
                saveas(hFig(i),[dest_folder,'Electro-JL-ACE-' prefix '-' num2str(j)],'epsc');
            else % F0mod plots
                plot(xaxis_tones(j)-50+H_dB_F0m(:,j)*50,freq_warped,'r')
                legend('seq','FFT')
                saveas(hFig(i),[dest_folder,'Electro-JL-F0m-' prefix '-' num2str(j)],'epsc');
            end
        end
    end
end

rmpath([directory 'PR_Stimuli-JL/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end