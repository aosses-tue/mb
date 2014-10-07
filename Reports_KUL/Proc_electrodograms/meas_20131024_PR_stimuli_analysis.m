function meas_20131024_PR_stimuli_analysis
% Programmed by Alejandro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'ConvertSequence_txt2nsb.m', ...
                    'ConvertSequence_nsb2NMT.m', ...
                    'rmsdb.m', ...
                    'Plot_sequence.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bSave = 1;
directory = '/home/alejandro/Documenten/Meas/Meas/Electrodograms/20131024-checking-RP-stimuli/';
addpath(directory)
addpath([directory 'PR_Stimuli-RP/'])
load([directory 'tg_p-q-structures_ACE_Simulink-var.mat']);

seqACE      = read_aseq([directory 'Sequence-104Hz-RP-ACE-highest.aseq']);
seqF0m      = read_aseq([directory 'Sequence-104Hz-RP-F0m-highest.aseq']);
seqACE_low  = read_aseq([directory 'Sequence-104Hz-RP-ACE-lowest.aseq']);
seqF0m_low  = read_aseq([directory 'Sequence-104Hz-RP-F0m-lowest.aseq']);

seqACEpatient = MapElectrodogram2PatientBoudaries(seqACE,p);
seqF0mpatient = MapElectrodogram2PatientBoudaries(seqF0m,p);
seqACEpatient_low = MapElectrodogram2PatientBoudaries(seqACE_low,p);
seqF0mpatient_low = MapElectrodogram2PatientBoudaries(seqF0m_low,p);

H_dB        = [];
H_dB_F0m    = [];
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_ACE_104_Hz.wav');
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_ACE_110_Hz.wav');
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_ACE_117_Hz.wav');
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_ACE_124_Hz.wav');
[freq_warped, H_dB(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_ACE_131_Hz.wav');

[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_F0m_104_Hz.wav');
[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_F0m_110_Hz.wav');
[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_F0m_117_Hz.wav');
[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_F0m_124_Hz.wav');
[freq_warped, H_dB_F0m(:,end+1)] = ConvertFFT2FAT('PR_Stimuli-RP/UW_LB_F0m_131_Hz.wav');

xaxis_tones     = [0 650 1380 2100 2830];
PlotSize        = [99 274 1811 425];
hFig    = [];
hAx     = [];
Plot_sequence(seqACEpatient)
title('ACE, balanced stimuli for RP + 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hACE_high = gca;
 
Plot_sequence(seqACEpatient_low)
title('ACE, balanced stimuli for RP - 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hACE_low = gca;

Plot_sequence(seqF0mpatient)
title('F0mod, balanced stimuli for RP + 5 dB')
hFig(end+1) = gcf;
hAx(end+1) = gca;
hF0m_high = gca;

Plot_sequence(seqF0mpatient_low)
title('F0mod, balanced stimuli for RP - 5 dB')
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
            if i <=2
                plot(xaxis_tones(j)-50+H_dB(:,j)*50,freq_warped,'r')
                legend('seq','FFT')
                saveas(hFig(i),[dest_folder,'Electro-RP-ACE-' prefix '-' num2str(j)],'epsc');
            else
                plot(xaxis_tones(j)-50+H_dB_F0m(:,j)*50,freq_warped,'r')
                legend('seq','FFT')
                saveas(hFig(i),[dest_folder,'Electro-RP-F0m-' prefix '-' num2str(j)],'epsc');
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end