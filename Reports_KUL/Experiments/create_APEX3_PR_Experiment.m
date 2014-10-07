function create_APEX3_PR_Experiment(ref_tone, outputfile, list_files, list_cal,bUseLB)
% function create_APEX3_PR_Experiment(ref_tone, outputfile, list_files, list_cal,bUseLB)
%
% ref_tone      - reference tone
% ref_octave    - number of octave
% outputfile    - file name for the created experiment
%
% For a version of this file where the L and R soundcard outputs were used,
% go back to a version committed before 19-08-2013. This is also valid for
% the file 'Template_PR.xml'
%
% % Example:
%       ref_tone    = 'Gsh'; 
%       ref_octave  = 3
%       outputfile  ='/home/alejandro/Documenten/Meas/Meas/Experiments/Tests_XML/Simple_demo_Tom_createdPR.xml'; 
%       create_APEX3_PR_Experiment(ref_tone, ref_tone, outputfile);
%
% All the wav files are mono files.
%
% Edit this function for changing the template and the location and names 
% of the audio files.
%
% Not yet:
%   Check generation of wav-files
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'note2num.m', ...
                    'num2note.m', ...
                    'setupPlotConf.m'};
if isunix
    UserName = 'alejandro';
else
    UserName = 'r0366612';
end

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files, UserName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isunix
    dirBase = ['C:\Documents and Settings\' UserName '\Desktop\Meas'];
else
    dirBase = ['/home/' UserName '/Documenten/Meas/Meas'];
end

template        =[dirBase delim 'Experiments' delim 'Tests_XML' delim 'Template_PR.xml'];
dirAudioFiles = '../../../Music/PR_Stimuli'; % Automate this respect to the Outputfile location

list_stim = {[list_files{2},'.wav'], ... % Automate audio file generation eliminating Frequency Number
             [list_files{3},'.wav'], ...
             [list_files{4},'.wav'], ...
             [list_files{5},'.wav']};

list_std = {[list_files{1},'.wav']};

list_cal = {[list_cal{1},'.wav']};

pr.trials       = ''; % structure containing what I want to write
pr.datablocks   = '';
pr.stimuli      = '';
%varParam.left   = '0'; % '0' is the value of the variable parameter
datablockcalid  = 'Noise'; % Noise filename
stimuluscalid   = 'calibrationstimulus';

%rove4cal        = varParam;
trial_id_number = 1:4; % 4 tones + reference = 5 tones

count                       = 0;
numRoving                   = 11; % Roving between +-5
step_dB                     = -1; % 0:step_dB:step_dB*numRoving (0:-3:-9) = 9 dB, to calibrate 0 dB @ 64 dB
pr.calibrationAmplitude     = 65;
pr.targetAmplitude          = 60;
numTrials                   = 60; % Ideal number of trials

pr.numPresentations         = numRoving;
num_ref_tone                = note2num(ref_tone.note);

for j = 1
    sufix_std_name      = ['_',ref_tone.note,'_',num2str(ref_tone.octave),'_',num2str(j)];
    standardid{j}       = sprintf('referencestimulus%s','');
    datablockstdid{j}   = sprintf('rawdata%s' , '_ref'); % 1 is for reference
    list_std_rov{j}     = list_std{1}; % just one reference
    
    for i=1:length(trial_id_number)
        num_tone        = mod(num_ref_tone + i, 12);
        if num_tone < num_ref_tone && num_tone ~= 0
            note_octave = ref_tone.octave + 1;
        else
            note_octave = ref_tone.octave;
        end

        tone = num2note(num_tone); % num2note edited for adding case 0 = case 12 = 'B'
        if strcmp( tone(end),'#')
            tone(end)=[];
            tone = [tone 'sh']; % from # to sh (sharp)
        end
        count = count + 1;
        trialid         = sprintf('trial%s', ['_',num2str(count)]);
        sufix_name       = ['_',tone    ,'_',num2str(note_octave),'_',num2str(j)];
        sufix_name_block = ['_',tone    ,'_',num2str(note_octave),'_',num2str(j)];
        stimulusid{count}   = sprintf('stimulus%s', num2str(i) ); % always 'rove 1' (roving being done by setting amp gain)
        datablockstimid{count} = sprintf('rawdata%s' , ['_test' num2str(i)] ); % 1 is for reference
        list_stim_rov{count} = list_stim{i};
        
        pr.trials       = [pr.trials a3trial(trialid, 'screen1', stimulusid{count},'',standardid{j}) ];
    end
end

disp(['Total created trials: ' num2str(count)])

pr.datablocks   = [pr.datablocks  a3datablocks([list_std_rov   list_stim_rov   list_cal],   dirAudioFiles, ...
                                               [datablockstdid datablockstimid datablockcalid])];

pr.stimuli      = [pr.stimuli     a3stimuli([datablockstdid datablockstimid datablockcalid], ...
                                            [standardid     stimulusid      stimuluscalid])]; %,[rove4std rove4stim rove4cal],'')];

output=readfile_replace(template, pr);
try
    fid=fopen(outputfile, 'w');
    fwrite(fid, output);
    fclose(fid);
catch
    mkdir([dirBase delim 'Experiments' delim 'Tests_XML' delim 'PR_new' delim]);
    disp(['Output directory: ' dirBase delim 'Experiments' delim 'Tests_XML' delim 'PR_new' delim])
    fid=fopen(outputfile, 'w');
    fwrite(fid, output);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end