function create_APEX3_LB_Experiment(ref_tone, tone, outputfile, list_files, list_cal)
% function create_APEX3_LB_Experiment(ref_tone, tone, outputfile, list_files, list_cal)
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
% Not yet:
%   Check generation of wav-files
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
    dirBase = ['/home/' UserName '/Documenten/Meas/Meas'];
end

template        =[dirBase delim 'Experiments/Tests_XML/Template_LB.xml'];
dirAudioFiles = '../../../Music/LB_Stimuli'; % Automate this respect to the Outputfile location

list_ref  = {[list_files{1},'.wav'] }; % Automate audio file generation eliminating Frequency Number
list_test = {[list_files{2},'.wav']};
list_cal  = {[list_cal{1},'.wav']};

pr.trials       = ''; % structure containing what I want to write
pr.datablocks   = '';
pr.stimuli      = '';
varParam.left   = '0'; % '0' is the value of the variable parameter
varParam.right  = '0';
datablockcalid  = 'Noise'; % Noise filename
stimuluscalid   = 'calibrationstimulus';

varGain         = varParam;
trial_id_number = 1;

count           = 0;

% To replace in Template:
pr.numPresentations = 1;
pr.nreversals       = 10;
pr.maxvalueac       = 50;
pr.step1            = 4;
pr.step2            = 2;
pr.file_ref         = list_ref{1};
pr.file_stim        = list_test{1};
pr.file_ref_cal     = list_cal{1};
pr.reference_level  = 0; %-10; % calibrated for 65 dB(A) then 65 dB(A) - 10 dB
pr.reversals_for_mean = 6;

if ref_tone.freq < 250
    ref_dBSPL       = 75; % at 131 Hz 75 dB SPL = 60 dB(A)
else
    ref_dBSPL       = 60;
end

try
    StartPoint      = round([NeuLoud(ref_tone.freq,ref_dBSPL)-NeuLoud(tone.freq,ref_dBSPL)]);
    display(['Inside ' mfilename])
    display(['Using starting point of ' num2str(StartPoint) ' respect to ' num2str(ref_dBSPL) ' dB SPL @ ' num2str(ref_tone.freq) ' Hz'])
    if ref_tone.freq == 131 && ref_dBSPL == 75
        display('75 dB SPL are approx. 60 dB(A)')
    end
catch
    StartPoint      = 0;
    display(['Inside ' mfilename])
    display('Function NeuLoud not found, we are going to start with 0 dB respect to the cal tone level')
    display('Press any button to continue...')
    pause()
end
    
deltaStartValue     = [StartPoint+10 StartPoint-10];

if deltaStartValue(1) > 5
    deltaStartValue(1) = 5; % referenced to pr.reference_level = -10
end

num_ref_tone    = note2num(ref_tone.note);

j = 1;

varParam.left   = 0; 
varParam.right  = varParam.left;

for k=1:length(deltaStartValue)
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
        stimulusid{count}   = sprintf('stimulus%s', sufix_name  ); % always 'rove 1' (roving being done by setting amp gain)
        datablockstimid{count} = sprintf('rawdata%s' , sufix_name_block ); % 1 is for reference
        varParam.left   = num2str(0); 
        varParam.right  = varParam.left;
        
        pr.StartValue   = deltaStartValue(k);
        try
            output      = readfile_replace(template, pr);
        catch
            error(['Error while generating XML experiments. Most probably the Template file was not found. For LB the template should be ' template])
        end
        
        if pr.StartValue >= 0
            output_gain = [outputfile(1:end-4) '-p' num2str(    pr.StartValue ) outputfile(end-3:end)];
        else
            output_gain = [outputfile(1:end-4) '-m' num2str(abs(pr.StartValue)) outputfile(end-3:end)];
        end
        
        try
            fid=fopen(output_gain, 'w');
            fwrite(fid, output);
            fclose(fid);
        catch
            mkdir([dirBase delim 'Experiments' delim 'Tests_XML' delim 'LB_new' delim]);
            fid=fopen(output_gain, 'w');
            fwrite(fid, output);
            fclose(fid);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end