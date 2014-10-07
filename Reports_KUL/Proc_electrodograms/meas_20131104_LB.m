function meas_20131104_LB

directory = '/home/alejandro/Documenten/Meas/Meas/Electrodograms/20131024-checking-RP-stimuli/PR_Stimuli-RP/'
addpath(directory)

audio_files = {'UW_LB_ACE_104_Hz.wav', 'UW_LB_F0m_104_Hz.wav'};

i = 1;

info.FileName = [directory, audio_files{i}];
info.VariableName = info.FileName;
info.bSave = 0;
[x, Fs] = wavread(info.FileName);
% p.CFG.Fs = 15659.375;
% % Map at 1800 pps
% p.STIM.CLevels = [170 169 168 170 172 172 173 174 174 174 174 174 175 176 176 176 175 174 173 174 174 172];
% p.STIM.TLevels = [117 115 113 111 109 110 112 117 122 120 117 120 124 129 134 129 124 118 111 114 117 103];

try
    p = evalin('base','p'); % Brings q-structure
catch
    info.pcl_filename = 'CS_RP_1800ppsF0.pcl'; 
    info.def_filename = 'nsb_SP15R1F0.def'; 
    ParMate(info.def_filename,info.pcl_filename);
end
        
x  = resample(x,round(p.CFG.Fs), Fs);  

[h, pF0m, ResF0m, pACE, info, vACE, vF0m] = Run_NMT(x, info, p);

% pImplant = Ensure_implant_params; % CIC3, Romain Peeters
pImplant = Current_level_to_microamps_proc;
uA_F0m = Current_level_to_microamps_proc(pImplant, vF0m.current_levels);

uA_F0mod = zeros(22,length(vF0m.current_levels/8));
for i = 1:22
    idx = find(vF0m.electrodes==i);
    uA_F0mod(i,1:length(idx)) = uA_F0m(idx);
end

RMS_F0m = sqrt( 1/(length(vF0m.current_levels)/8)*sum(uA_F0mod' .* uA_F0mod'));

uA_ACE = Current_level_to_microamps_proc(pImplant, vACE.current_levels);

uA_ACEtot = zeros(22,length(vACE.current_levels/8));
for i = 1:22
    idx = find(vACE.electrodes==i);
    uA_ACEtot(i,1:length(idx)) = uA_ACE(idx);
end

RMS_ACE = sqrt( 1/(length(vACE.current_levels)/8)*sum(uA_ACEtot' .* uA_ACEtot'));

DeltaRMS = RMS_ACE - RMS_F0m;

Delta_CU = Microamps2current(DeltaRMS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end