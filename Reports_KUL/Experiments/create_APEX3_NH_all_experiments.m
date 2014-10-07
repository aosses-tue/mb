function create_APEX3_NH_all_experiments(bPart, bUseLB, diroutput)
% function create_APEX3_NH_all_experiments(bPart, bUseLB, diroutput)
%
% Referential nomenclature for apex result file names:
%
%   Experiment files: \textsf{XX\_Ref\_YYY\_ZZ.xml}\vspace{0pt}
%   Result files: \textsf{XX\_Ref\_YYY\_ZZ\_AAA-SS-N.apr}\vspace{12pt}
%
% Wav-file names for PR experiment:     'UW_104_Hz'         % without LB
%                                       'UW_LB_ACE_104_Hz'  % with LB for ACE
% Example 1:
%   bPart = 0; bLB = 1; create_APEX3_NH_all_experiments(bPart,bLB);
%
% Example 2: Only LB
%   bPart = 3; bLB = 0; create_APEX3_NH_all_experiments(bPart,bLB);
%
% Example 3: Generating PR and MCI with LB:
%   bPart = 0; bLB = 1; create_APEX3_NH_all_experiments(bPart,bLB);
%   
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'cos_ramp.m',...
                    'eq_rms', ...
                    'get_tones_for_pitch_ranking.m',...
                    'mweightingfilter.m',...
                    'NeuLoud.m', ... % Defines starting point for loudness balance
                    'note2freq.m',...
                    'UW_CAMP_STIMULI_AMP.mat'}; 
                
% create_uw_camp_stimuli_new
% Not committed yet:
%   create_from_raw_stimuli  

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files,'alejandro');

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

diroutput = [dirBase delim 'Experiments' delim 'Tests_XML' delim];

if bUseLB == 1
    times_LB_procedure = 2;
else
    times_LB_procedure = 1;
end

bGenerateTones  = 1; % Set to 0 when the wav-files are already created

instr           = {'UW'};   % Just one instrument according to meeting with TF on 31/05/2013

display('Created experiments:')

[idle, tones]   = get_tones_for_pitch_ranking(0);

if bUseLB == 0
    [wavfilenames  , calfilenames  , DirAudioFilesPR] = create_PR_NH_stimuli(idle, tones, instr, 0); % Just get the names for Wav-files
    [wavfilenamesLB, calfilenamesLB, DirAudioFilesLB] = create_LB_NH_stimuli(idle, tones, instr, 1); % 1 - bLB = bGenerateTones for LB. 

    wavfilenamesLB      = { wavfilenames{1,1} wavfilenames{2,3} wavfilenames{4,1} wavfilenames{5,3} };
    basefilenameLB      = { wavfilenames{2,1} };
    basetonesLB         = { idle{2,1} }; % Always the tone of 131 Hz
    comparisontonesLB   = { idle{1,1} idle{2,3} idle{4,1} idle{5,3} };

    if bPart == 0 || bPart == 3

        diroutput_LB    = [diroutput 'LB_new' delim];
        %% Loudness Balance
        XX  = 'LB';

        % Comparison tones = 104, 147, 208, 294 Hz
        % Reference tone = 131 Hz

        for j = 1:length(instr)
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
end



[idle, tones]   = get_tones_for_pitch_ranking(0);
[wavfilenames  , calfilenames  , DirAudioFilesPR, LB_Gains] = create_PR_NH_stimuli(idle, tones, instr, bGenerateTones, bUseLB); % used in PR but also in MC

if bPart == 0 || bPart == 1
    
    diroutput_PR    = [diroutput 'PR_new' delim];
    %% Pitch Ranking
    Strategy = [];    
    XX  = 'PR';
    
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
    Strategy = [];
    
    for m = 1:times_LB_procedure
    
        if m == 1 
            if bUseLB == 1
                Strategy = '_ACE';
            end
        else
            Strategy = '_F0m';
        end
        
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
        conf.ints       = {'1','2','3','4'};
        conf.Strategy   = Strategy;
        
        wavfilenames    = create_MC_NH_stimuli(ref_tone, conf, bGenerateTones, bUseLB, LB_Gains_Strategy); % Reusing Cal-files
        [a, b]          = size( wavfilenames );
        CantIntervals   = a / length(instr); % CantIntervals is similar to the variable Num_series


        for j = 1:length(instr)

            ZZ          = instr{j}; % instrument
            for i = 1:CantIntervals

                wav_strings = wavfilenames(i+(j-1)*CantIntervals,:); 
                calfilename = calfilenames(j);

                filename    = [XX,'_Ref_', num2str( round(note2freq(ref_tone)) ),'_Int_',num2str(i),'-',ZZ,Strategy,'.apx'];
                outputfile  = [diroutput_MC, filename]; % Each interval, each instrument
                create_APEX3_MC_Experiment(ref_tone, outputfile, wav_strings, calfilename, bUseLB, Strategy);

                display(['file name: ',filename]);
            end

        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end