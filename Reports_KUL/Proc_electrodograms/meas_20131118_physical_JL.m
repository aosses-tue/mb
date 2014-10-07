function meas_20131118_physical_JL
% function meas_20131118_physical_JL
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cal-JL-xPC-ACE-v097.out
% /media/Elements/orl-wrk-0089/Documenten/Meas/Meas/Electrodograms/20131118-Physical-JL-stimuli/tg_p-q-structures_Simulink-var.mat

% Programmed by Alejandro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'ConvertSequence_txt2nsb.m', ...
                    'ConvertSequence_nsb2NMT.m', ...
                    'Get_average_current_levels.m', ...
                    'rmsdb.m', ...
                    'Plot_sequence.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bSave = 1;
directory = '/media/Elements/orl-wrk-0089/Documenten/Meas/Meas/Electrodograms/20131118-Physical-JL-stimuli/';
addpath(directory)
load([directory 'tg_p-q-structures_Simulink-var.mat']); % Subject map

try
    XLim = [1000 1050];
    seq = read_aseq([directory 'cal-JL-xPC-ACE-v097.aseq']); % Sequence recorded with the xPC system
    Plot_sequence_patient(seq,p);
    title('xPC ACE, Pure tone at 1 kHz, map: JL')
    xlim(XLim)
    
    idx_min = Get_num_pulses_up_to(seq,XLim(1)/1000);
    idx_max = Get_num_pulses_up_to(seq,XLim(2)/1000);
    Avg = Get_average_current_levels(seq, idx_min, idx_max,p.STIM.NumMaxima);
    
catch
    seq = ConvertSequence_txt2nsb('cal-JL-xPC-ACE-v097.out');
    [seq seqPatient] = ConvertSequence_nsb2NMT(seq,p);
    % save_aseq(seq)
    Plot_sequence_patient(seq,p);
    title(['Map: JL, Volume = ' num2str(p.STIM.Volume) '; Strategy: xPC ACE'])
end

seqACE = read_aseq([directory 'PureTone-1kHz-57dB-ACE-Simulink-p08.aseq']);

Plot_sequence_patient(seqACE,p);
title('Map: JL; Strategy: Sim ACE')

XLim = [1000 1050];
xlim(XLim)
idx_min_ACE = Get_num_pulses_up_to(seqACE,XLim(1)/1000);
idx_max_ACE = Get_num_pulses_up_to(seqACE,XLim(2)/1000);
AvgACE = Get_average_current_levels(seqACE, idx_min_ACE, idx_max_ACE,p.STIM.NumMaxima);

disp('Elec.Num. Avg-xPC Avg-Sim')
disp(num2str([(1:22)' Avg AvgACE]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end