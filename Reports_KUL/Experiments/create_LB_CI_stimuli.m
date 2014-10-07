function [wavfilename calfilename DirAudioFiles] = create_LB_CI_stimuli(idle, tones, instr, bGenerateTones, bA_weighting)
% function [wavfilename calfilename DirAudioFiles] = create_LB_CI_stimuli(idle, tones)
%
% As an example of the use of this function, run the script
% create_APEX3_CI_all_experiments.
%
% The roving is being introduced in APEX3
%
% Dependencies:
%   get_tones_for_pitch_ranking
%   create_uw_camp_stimuli
%   rmsdb
%
% Programmed by Alejandro but adapted from Matthias' script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    bGenerateTones = 1;
end

if nargin < 5
    bA_weighting    = 0;
end

wavfilenames    = [];

cur_dir_old     = cd;

if ~isunix
    dirBase = 'C:\Documents and Settings\r0366612\Desktop\Meas';
else
    dirBase = '/home/alejandro/Documenten/Meas/Meas';
end

cur_dir_new     = [dirBase delim 'Music' delim 'LB_Stimuli_new' delim];
try
    cd(cur_dir_new);
catch
    mkdir(cur_dir_new)
    disp(['Inside ' mfilename ' : ' cur_dir_new ' created. Wav-files will be generated here']);
    cd(cur_dir_new);
end

DirAudioFiles = cur_dir_new;

if nargin == 0
    [idle, tones] = get_tones_for_pitch_ranking(0);
end

[a,b]       = size(idle); 
Num_tones_per_serie = b;
Num_series  = a;

fs          = 44100;
stim_len    = 0.5; % in seconds
cal_len     = 3; % in seconds, for calib
ramp_len    = 25; % in mili seconds
ramp_type   = '2'; % cos rising and falling ramp
off         = 0.001;
cont        = 0;
dBFS        = -20; % -20 dBFS
freq_cal    = 131;

if bA_weighting
    display('A-Weighting being applied...')
    A_Curve = mweightingfilter('A', fs, 0);
end

%% Calibration Tone:

for k=1:length(instr)
    
    for i = 5 % i = 5 corresponds to calibration tone
        
        if bA_weighting == 0
            calfilenames    = {[ instr{k} '_' num2str(freq_cal) '_Hz-cal']};
        end
        
        if bGenerateTones
            tmp         = strread(tones{1,i}, '%s', 'delimiter', '_');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % After SVN revision 410 only UW filters being used
         	multisine   = create_uw_camp_stimuli(           freq_cal, fs, cal_len+2*ramp_len*1e-3, ramp_type, ramp_len ); % 0 ms
            idxEliminateFadeInOut = ceil(ramp_len*1e-3*fs);
            multisine   = multisine(idxEliminateFadeInOut : end-idxEliminateFadeInOut);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            crms        = rmsdb(multisine); % RMS value of the stimulus just after created
            if bA_weighting
                amultisine = filter(A_Curve.numer, A_Curve.denom, multisine); % A-weighted Sine tone
                arms    = rmsdb(amultisine); % RMS value of the A-weighted stimulus
            else
                arms    = crms; % Thus crms - arms = 0
            end
            delta_A     = crms - arms;
            multisine   = eq_rms(multisine, dBFS+delta_A);
            wavwrite( multisine, fs, calfilenames{1});
            cont     	= cont + 1; % counting created WAV-files
        end
        
    end

    for i = 1:length( tones )

        tmp         = strread(tones{1,i}, '%s', 'delimiter', '_');
        tone.note   = tmp{1};
        tone.octave = str2num(tmp{2});
        if bA_weighting == 0
            wavfilenames{1,i} = [ instr{k},'_', tones{2, i} '_Hz'];
        end
        
        if bGenerateTones
            if strcmp( instr(k),'UW' )
                multisine   = create_uw_camp_stimuli(         note2freq( tone ), fs, stim_len, ramp_type, ramp_len );
            end
            crms        = rmsdb(multisine);
            if bA_weighting == 0
                arms    = crms;
            end
            delta_A     = crms-arms;
            multisine   = eq_rms(multisine, dBFS+delta_A);
            wavwrite( multisine, fs, wavfilenames{1,i});
            cont        = cont + 1; % counting created WAV-files
        end

    end

    for i = 1:4:length( tones )-1 
        aux_cont = floor(i/4)+1 + (k-1)*Num_series; 
        disp([mfilename '.m: reordering wav-filenames'])
        disp(['Row ' int2str(aux_cont) ' wav-file names from/to idx ' num2str(i) ' / ' num2str(i+Num_tones_per_serie-1)])
        wavfilename(aux_cont,1:Num_tones_per_serie) = wavfilenames(1,i:i+Num_tones_per_serie-1);
        calfilename(k,1) = calfilenames;
    end
end

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
if cont ~= 0 % Then audio files were already generated
    display([mfilename ': Wav files created into directory:'])
    display(cd)
    display([num2str(cont),' WAV-files',' at Fs=',num2str(fs),' Hz. Duration of_',num2str(stim_len),' seconds'])
    display('Please move the generated WAV-files to the corresponding folder')
else
    display([mfilename ': Wav files were not created, but you got the names of them, for creating the APEX Experiments'])
end
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd(cur_dir_old)
