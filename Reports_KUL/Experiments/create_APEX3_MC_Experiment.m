function create_APEX3_MC_Experiment(ref_tone, outputfile, list_files, list_cal, bUseLB, Strategy)
% function create_APEX3_MC_Experiment(ref_tone, outputfile, list_files, list_cal, bUseLB, Strategy)
%
% ref_tone      - reference tone
% ref_octave    - number of octave
% outputfile    - file name for the created experiment
%
% For instance:
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
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'note2num.m', ...
                    'num2note.m'};
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
    dirBase = ['/home/' UserName '/Documenten/Meas/Meas/'];
end

template        = [dirBase delim 'Experiments' delim 'Tests_XML' delim 'Template_MC.xml'];
dirAudioFiles   = '../../../Music/MC_Stimuli'; 

for i = 1:length(list_files)
    list_stim{i} = [list_files{i},'.wav'];
end
             
list_cal = {[list_cal{1},'.wav']};

pr.picturesFolder = '../../Pictures';
pr.trials       = ''; % structure containing what I want to write
pr.datablocks   = '';
pr.stimuli      = '';

if bUseLB
    pr.buttonLB = ['_LB' Strategy];
end

varParam.left   = '0'; % '0' is the value of the variable parameter
datablockcalid  = 'Noise'; % Noise filename
stimuluscalid   = 'calibrationstimulus';

casacal         = varParam;
trial_id_number = 1:length( list_files ); % 4 tones + reference = 5 tones

count       = 0;
numRoving   = 5; % 3
step_dB     = -2; % -2
pr.calibrationAmplitude     = 60 + abs(step_dB)*(numRoving-1)/2;
pr.targetAmplitude          = pr.calibrationAmplitude;

numTrials   = 60;
 
pr.numPresentations = floor(numTrials/length(trial_id_number)/numRoving);

num_ref_tone    = note2num(ref_tone.note);

for j = 1:numRoving
        
    for i=1:length(trial_id_number)
        count = count + 1;
        trialid         = sprintf('trial%s', ['_',num2str(count)]);
        if strcmp( list_stim{i}(3:5), '-A-' )
            seq     = list_stim{i}(6:end);
        else
            seq     = list_stim{i}(4:end);
        end
        seq     = seq(1:end-11); % deleting '_104_Hz.wav'
        answer  = ['button_',seq(1:end-1)];
        
        sufix_name       = ['_',seq,'_',num2str(j)];
        sufix_name_block = ['_',seq,'_',num2str(j)];
        stimulusid{count}   = sprintf('stimulus%s', sufix_name  ); % always 'rove 1' (roving being done by setting amp gain)
        datablockstimid{count} = sprintf('rawdata%s' , sufix_name_block ); % 1 is for reference
        varParam.left   = num2str((j-1)*step_dB); 
        casastim{count} = varParam;
        list_stim_rov{count} = list_stim{i};
        
        pr.trials       = [pr.trials a3trial(trialid, 'screen1', stimulusid{count}, answer) ];
    end
end

disp(['Total created trials: ' num2str(count*pr.numPresentations)])

pr.datablocks   = [pr.datablocks    a3datablocks(   [list_stim_rov   list_cal],   dirAudioFiles, ...
                                                    [datablockstimid datablockcalid])];

pr.stimuli      = [pr.stimuli       a3stimuli(      [datablockstimid datablockcalid], ...
                                                    [stimulusid      stimuluscalid],[casastim casacal],'')];

output=readfile_replace(template, pr);

try
    fid=fopen(outputfile, 'w');
    fwrite(fid, output);
    fclose(fid);
catch
    mkdir([dirBase delim 'Experiments' delim 'Tests_XML' delim 'MC_new' delim]);
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
