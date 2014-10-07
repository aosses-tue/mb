function create_APEX3_CI_all_experiments(bPart, bUseLB, diroutput)
% function create_APEX3_CI_all_experiments(bPart, bUseLB, diroutput)
%
% Referential nomenclature for apex result file names:
%
%   Experiment files: \textsf{XX\_Ref\_YYY\_ZZ.xml}\vspace{0pt}
%   Result files: \textsf{XX\_Ref\_YYY\_ZZ\_AAA-SS-N.apr}\vspace{12pt}
%
% Wav-file names for PR experiment:     'UW_104_Hz'         % without LB
%                                       'UW_LB_ACE_104_Hz'  % with LB for ACE
%
% Scripts where the experiment files are actually being generated:
%   LB experiments:  create_APEX3_LB_Experiment.m
%   PR experiments:  create_APEX3_PR_Experiment.m
%   MCI experiments: create_APEX3_MC_Experiment.m
%
% Example 1:
%   bPart = 0; bLB = 1; create_APEX3_CI_all_experiments(bPart,bLB);
%
% Example 2: Only LB
%   bPart = 3; bLB = 0; create_APEX3_CI_all_experiments(bPart,bLB);
%
% Example 3: Generating PR and MCI with LB:
%   bPart = 0; bLB = 1; create_APEX3_CI_all_experiments(bPart,bLB);
%   
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'cos_ramp.m',...
                    'eq_rms', ...
                    'get_tones_for_pitch_ranking.m',...
                    'NeuLoud.m', ... % Defines starting point for loudness balance
                    'note2freq.m',...
                    'UW_CAMP_STIMULI_AMP.mat'}; 
                
% create_uw_camp_stimuli_new
% Not committed yet:
%   create_from_raw_stimuli  

if isunix ==  1
    UserName = 'alejandro';
else
    UserName = 'r0366612' ;
end

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files,UserName);

if isunix
    addpath('/home/alejandro/Documenten/MM/src/Matlab/Classes/')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bPart = 0;
end

if nargin < 2
    bUseLB = 0;
end

if nargin < 3
    if ~isunix
        dirBase   = 'C:\Documents and Settings\r0366612\Desktop\Meas';
    else
        dirBase   = '/home/alejandro/Documenten/Meas/Meas';
    end
end

subject_type = 'CI';

diroutput = [dirBase delim 'Experiments' delim 'Tests_XML' delim];

if bUseLB == 1
    times_LB_procedure = 2;
else
    times_LB_procedure = 1;
end

bGenerateTones  = 1; % Set to 0 when the wav-files are already created

instr           = {'UW'};   % Just one instrument according to meeting with TF on 31/05/2013

display('Created experiments:')

[idle, tones]   = get_tones_for_pitch_ranking(2);

% Default values:
[wavfilenames  , calfilenames  , DirAudioFilesPR, LB_Gains]     = create_PR_CI_stimuli(idle, tones, instr, 0); % Just get the names for Wav-files

if bUseLB == 0
    [wavfilenamesLB, calfilenamesLB, DirAudioFilesLB]           = create_LB_CI_stimuli(idle, tones, instr, 1); % 1 - bLB = bGenerateTones for LB. 

    wavfilenamesLB      = { wavfilenames{1,1} wavfilenames{2,3} wavfilenames{4,1} ...
                            wavfilenames{5,3} wavfilenames{7,4} wavfilenames{9,5} }; % 'UW_294_Hz' 'UW_494_Hz' 'UW_831_Hz', respectively
    basefilenameLB      = { wavfilenames{2,1} };
    basetonesLB         = { idle{2,1} }; % Always the tone of 131 Hz
    comparisontonesLB   = { idle{1,1} idle{2,3} idle{4,1} idle{5,3} idle{7,4} idle{9,5} };

    if bPart == 0 || bPart == 3

        diroutput_LB    = [diroutput 'LB_new' delim];
        disp([mfilename '.m: the experiments will be generated at '' ' diroutput_LB ''''])
        %% Loudness Balance
        XX  = 'LB';

        % Comparison tones = 104, 147, 208, 294, 494, 831 Hz
        %                                        B4 , G#4
        % Reference tone = 131 Hz

        j           = 1; % Just one instrument
        ZZ          = instr{j}; % instrument
        Num_series  = length(basetonesLB);

        for i=1:length(comparisontonesLB)

            wav_strings = [basefilenameLB(1,1) wavfilenamesLB(1,i)]; 
            calfilename = calfilenamesLB(j);

            tones       = comparisontonesLB{1, i};
            tone.note   = tones(1:end-2);
            tone.octave = str2num(tones(end));
            tone.freq   = round( note2freq(tone) );

            tones       = basetonesLB{1,1}; % Always reference
            toneRef.note   = tones(1:end-2);
            toneRef.octave = str2num(tones(end));
            toneRef.freq   = round( note2freq(toneRef) );

            filename    = [XX,'_Com_', num2str(tone.freq),'_',ZZ,'.apx'];
            outputfile  = [diroutput_LB, filename]; 
            if bUseLB == 0
                create_APEX3_LB_Experiment(toneRef, tone, outputfile, wav_strings, calfilename);
            end

            display(['file name: ',filename]);
        end
    end
end

if bPart == 0 || bPart == 1
    
    diroutput_PR    = [diroutput 'PR_new' delim];
    disp([mfilename '.m: the experiments will be generated at '' ' diroutput_PR ''''])
    %% Pitch Ranking
    Strategy = [];    
    XX  = 'PR';
    
    if strcmp(subject_type,'CI')
        [idle, tones]   = get_tones_for_pitch_ranking(2);
    else
        [idle, tones]   = get_tones_for_pitch_ranking(0);
    end
        
    [wavfilenames  , calfilenames  , DirAudioFilesPR, LB_Gains, test_regs] = create_PR_CI_stimuli(idle, tones, instr, bGenerateTones, bUseLB, subject_type);

    basetones       = idle(: , 1); % First column
    Num_series      = length(basetones);
    
	for m = 1:times_LB_procedure
        
        if m == 1 
            if bUseLB == 1
                Strategy = '_ACE';
            end
        else
            Strategy = '_F0m';
        end
        
        wavfilenamesPR  = wavfilenames((1:Num_series)+(m-1)*Num_series,:);
        calfilenamePR   = calfilenames(m);
        
        for j = 1:length(instr)
            ZZ          = instr{j}; % instrument

            for i=1:length(basetones)

                wav_strings = wavfilenamesPR(i+(j-1)*Num_series,:); 
                calfilename = calfilenamePR;

                tones       = basetones{i};
                tone.note   = tones(1:end-2);
                tone.octave = str2num(tones(end));
                tone.freq   = round( note2freq(tone) );
                if bUseLB == 0
                    filename    = [XX,'_Ref_', num2str(tone.freq), '_', ZZ, Strategy, '.apx'];
                else
                    filename    = [XX,'_Ref_', num2str(tone.freq), '_', ZZ, Strategy, '-LB.apx'];
                end
                outputfile  = [diroutput_PR, filename]; 
                create_APEX3_PR_Experiment(tone, outputfile, wav_strings, calfilename,bUseLB);

                display(['file name: ',filename]);
            end
        end
    end
    
end

%% MCI
XX  = 'MC';

% We copy the same calibration tone used in the PR Experiment
if bPart == 0;
    
    for m = 1:times_LB_procedure
        srcCalFileName  = [DirAudioFilesPR char(calfilenames(m)) '.wav'];
        DirAudioFilesMC = dir_change_level(DirAudioFilesPR, 1, 'MC_Stimuli');
        destCalFilename = [DirAudioFilesMC delim char(calfilenames(m)) '.wav'];

        copyfile(srcCalFileName, destCalFilename);
    end

end

if bPart == 0 || bPart == 2
    
    diroutput_MC    = [diroutput 'MC_new' delim];
    disp([mfilename '.m: the experiments will be generated at '' ' diroutput_MC ''''])
    
    Strategy = [];
    
    for m = 1:times_LB_procedure
    
        if m == 1 
            if bUseLB == 1
                Strategy = '_ACE';
            end
        else
            Strategy = '_F0m';
        end
        % size(LB_Gains) = 2 x 42; freqs from 104 Hz to 330 Hz
        if bUseLB
            if m == 1
                LB_Gains_Strategy = LB_Gains(:,1:end/2);
            else
                LB_Gains_Strategy = LB_Gains(:,end/2+1:end);
            end
        else
            
            LB_Gains_Strategy = LB_Gains; % Always null-values
            
        end
        
        ref_tone.note   = 'Gsh';
        ref_tone.octave = 2;
      % conf.ints       = {'1','2','3','4'}; % This is for NH
        conf.ints       = {'1','3','6','9'}; % This is for CI
        conf.Strategy   = Strategy;
        
        wavfilenames    = create_MC_NH_stimuli(ref_tone, conf, bGenerateTones, bUseLB, LB_Gains_Strategy); % Reusing Cal-files
        [a, b]          = size( wavfilenames );
        CantIntervals   = a / length(instr); % CantIntervals is similar to the variable Num_series


        for j = 1:length(instr)

            ZZ          = instr{j}; % instrument
            for i = 1:CantIntervals

                wav_strings = wavfilenames(i+(j-1)*CantIntervals,:); 
                calfilename = calfilenames(j);

                filename    = [XX,'_Ref_', num2str( round(note2freq(ref_tone)) ),'_Int_',conf.ints{i} ,'-',ZZ,Strategy,'.apx'];
                outputfile  = [diroutput_MC, filename]; % Each interval, each instrument
                create_APEX3_MC_Experiment(ref_tone, outputfile, wav_strings, calfilename, bUseLB, Strategy);

                display(['file name: ',filename]);
            end

        end
    end
    
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([mfilename '.m: script ending. All the experiments were generated'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end