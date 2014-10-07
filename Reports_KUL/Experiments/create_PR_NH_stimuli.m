function [wavfilename calfilename DirAudioFiles, LB_Gain] = create_PR_NH_stimuli(idle, tones, instr, bGenerateTones, bUseLB, bA_weighting)
% function [wavfilename calfilename DirAudioFiles, LB_Gain] = create_PR_NH_stimuli(idle, tones, instr, bGenerateTones, bLB, bA_weighting)
%
% As an example of the use of this function, run the script
% create_APEX3_NH_all_experiments.
%
%   bUseLB = 1; means that Loudness Balancing is going to be applied
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

if nargin < 4
    bGenerateTones = 1;
end

if nargin < 5
    bUseLB         = 0;
end

if nargin < 6
    bA_weighting = 0;
end

bAnswer = 0;

if ~isunix
    dirBase = 'C:\Documents and Settings\r0366612\Desktop\Meas';
else
    dirBase = '/home/alejandro/Documenten/Meas/Meas';
end

dirResults = [dirBase delim 'Experiments' delim 'Results_XML' delim];

if bUseLB == 1
    [aceData, f0mData, subjects, stLB] = figLBScores_xPC(dirResults,{'UW'});
    nSubjects = size(aceData,1);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display([int2str(nSubjects), ' subjects found. Please make a choice:'])
    for i = 1:nSubjects
        display(['    Type ',int2str(i), ' to load the data for ', subjects{i, 1}]);
    end
    disp(' ')
    bAnswer = input('Please make your choice. If you do not want to use any of the stored data, type 0: ');
    
    if bAnswer == 0
        stLB = getLBValues;
        LB_values = [stLB.LB_ACE; stLB.LB_F0m];
        display('Loading fixed values...');
    else
        LB_values = [stLB.LB_ACE; stLB.LB_F0m]; % 1 row for each strategy
    end
    times_LB_procedure = 2; % ACE and F0mod
else
    times_LB_procedure = 1;
end

if bAnswer == 0 && bUseLB == 0; % Then no LB
    [aceData, f0mData, subjects, stLB] = figLBScores_xPC(dirResults,{'UW'});
    LB_values = [0*stLB.LB_ACE(1,:); 0*stLB.LB_F0m(1,:)];
end
% Remember that if LB is desired in a manual way...then bUseLB = 1; bAnswer
% = 0

LB_Gain         = [];
wavfilenames    = [];

cur_dir_old     = cd;
cur_dir_new     = [dirBase delim 'Music' delim 'PR_Stimuli' delim];
try
    cd(cur_dir_new);
catch
    mkdir(cur_dir_new)
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

Add2Name_if_Corrections = [];
Strategy = [];

if bA_weighting == 1
    Add2Name_if_Corrections = [Add2Name_if_Corrections '-A'];
end

if bUseLB == 1
    Add2Name_if_Corrections = [Add2Name_if_Corrections '_LB'];
end

Add2Name_if_Corrections = [Add2Name_if_Corrections '_']; % Default name for this variable: '_'

%% Calibration Tone:

for j = 1:times_LB_procedure
    
    if j == 1 && bUseLB == 1
        Strategy = 'ACE_';
    end
    if j == 2
        Strategy = 'F0m_';
    end
    
    for k=1:length(instr)

        for i = 5 % i = 5 corresponds to calibration tone

            calfilenames    = {[ instr{k} Add2Name_if_Corrections Strategy num2str(freq_cal) '_Hz-cal']};

            if bGenerateTones
                tmp         = strread(tones{1,i}, '%s', 'delimiter', '_');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % After SVN revision 410 only UW filters being used
                multisine = create_uw_camp_stimuli(           freq_cal, fs, cal_len+2*ramp_len*1e-3, ramp_type, ramp_len );
                idxEliminateFadeInOut = ceil(ramp_len*1e-3*fs);
                multisine   = multisine(idxEliminateFadeInOut : end-idxEliminateFadeInOut);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                crms        = rmsdb(multisine); % RMS value of the stimulus just after created
                if bA_weighting
                    amultisine = filter(A_Curve.numer, A_Curve.denom, multisine); % A-weighted Sine tone
                    arms    = rmsdb(amultisine); % RMS value of the A-weighted stimulus
                else
                    arms    = crms; % Thus crms - arms = 0
                end
                if bUseLB == 0 % i.e. no loudness balance
                    delta_A     = crms - arms; % Traditional calculation
                else
                    bFreqMatch  = 0;
                    m           = 1;
                    while bFreqMatch == 0
                        if stLB.F(m) == freq_cal;
                            delta_A = round(LB_values(j,m)); 
                            bFreqMatch = 1;
                        end
                        m = m+1;
                    end
                end
                multisine   = eq_rms(multisine, dBFS+delta_A);
                wavwrite( multisine, fs, calfilenames{1});
                cont        = cont + 1; % counting created WAV-files
            end

        end

        for i = 1:length( tones )

            tmp         = strread(tones{1,i}, '%s', 'delimiter', '_');
            tone.note   = tmp{1};
            tone.octave = str2num(tmp{2});
            freq_tone   = str2num(tones{2, i});
            wavfilenames{1,i} = [ instr{k},Add2Name_if_Corrections, Strategy, tones{2, i} '_Hz'];

            if bGenerateTones
                if strcmp( instr(k),'UW' )
                    multisine   = create_uw_camp_stimuli( freq_tone, fs, stim_len, ramp_type, ramp_len );
                end
                crms        = rmsdb(multisine);
                if bA_weighting
                    amultisine = filter(A_Curve.numer, A_Curve.denom,multisine);
                    arms    = rmsdb(amultisine);
                else
                    arms    = crms;
                end
                if bUseLB == 0 % i.e. no loudness balance
                    delta_A     = crms - arms; % Traditional calculation
                    bFreqMatch  = 0;
                    m           = 1;
                    while bFreqMatch == 0
                        if stLB.F(m) == freq_tone;
                            LB_Gain = [ LB_Gain [delta_A;freq_tone] ]; % Using null value
                            bFreqMatch = 1;
                        end
                        m = m+1;
                    end
                else
                    bFreqMatch  = 0;
                    m           = 1;
                    while bFreqMatch == 0
                        if m <= length(stLB.F)
                            if stLB.F(m) == freq_tone;
                                delta_A = LB_values(j,m); 
                                LB_Gain = [ LB_Gain [delta_A;freq_tone] ];
                                bFreqMatch = 1;
                            end
                            m = m+1;
                        else
                            bFreqMatch = 1;
                        end
                    end
                end
                multisine   = eq_rms(multisine, dBFS+delta_A);
                wavwrite( multisine, fs, wavfilenames{1,i});
                cont        = cont + 1; % counting created WAV-files   
            end

        end

        for i = 1:4:length( tones )-1 
            aux_cont = floor(i/4)+1 + (k-1)*Num_series; 
            wavfilename(aux_cont+(j-1)*Num_series,1:Num_tones_per_serie) = wavfilenames(1,i:i+Num_tones_per_serie-1);
            calfilename(j,1) = calfilenames;
        end
    end
end

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
if cont ~= 0 % Then audio files were already generated
    display([mfilename ': Wav files created into directory:'])
    display(cd)
    display(strcat(num2str(cont),' WAV-files','_at Fs=',num2str(fs),' Hz. Duration of_',num2str(stim_len),' seconds'))
    display('Please move the generated WAV-files to the corresponding folder')
else
    display([mfilename ': Wav files were not created, but you got the names of them, for creating the APEX Experiments'])
end
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

cd(cur_dir_old)
