function plot_electrodogram_with_fft
% function plot_electrodogram_with_fft
% 
% Adapted from: meas_20131025_PR_stimuli_analysis
%
% Programmed by Alejandro Osses
% Created in     : 21/05/2014
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

seqF0m      = read_aseq([directory 'Sequence-104Hz-JL-F0m-highest.aseq']);

seqF0mpatient = MapElectrodogram2PatientBoudaries(seqF0m,p);

H_dB_F0m    = [];
[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_104_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_110_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_117_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_124_Hz.wav');
[~          , H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-JL/UW_LB_F0m_131_Hz.wav');

xaxis_tones = [0 650 1380 2100 2830]; % Time in mili seconds where each tone starts

hFig    = [];
hAx     = [];
Plot_sequence(seqF0mpatient)
title('F0mod, balanced stimuli for JL + 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hF0m_high = gca;

PlotSize        = [99 274 1811 425];
for i = 1:length(hFig)
    set(hFig(i),'Position', PlotSize)
end


% hoofd_folder    = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/';
% dest_folder     = [hoofd_folder 'ci_pooled/Figures/'];

ylim([0 20])

for j = 1:length(xaxis_tones) % Set x min/max according to electrodogram position
    
    for i = 1:length(hFig) % 
        if bSave
            figure(hFig(i));
            xlim([xaxis_tones(j)-50 xaxis_tones(j)+600])
            rectangle('Position',[xaxis_tones(j)-50 0.01 50 20-0.05],'FaceColor','w'), hold on
            % Plots a white rectangle
            
            plot(xaxis_tones(j)-50+H_dB_F0m(:,j)*50,freq_warped,'r')
            legend('seq','FFT')
            %saveas(hFig(i),[dest_folder,'Electro-JL-F0m-' prefix '-' num2str(j)],'epsc');
        end
    end
end

rmpath([directory 'PR_Stimuli-JL/'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end