function experiment_training_20131119
% function experiment_training_20131119
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

location_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Tests_XML/';
directory       = 'Gen-speech-training/';
List_of_files   = {'oz_speech_material.m'}; 

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath([location_folder directory 'matlab-20131118/']);

p.directory = '/home/alejandro/Documenten/Meas/Meas/Experiments/Tests_XML/';

% p.speechmaterial = 'BKB';
p.speechmaterial    = 'VlMatrix';
ac_speech_training(p);

p.speechmaterial    = 'LISTf';
ac_speech_training(p);

p.speechmaterial    = 'LISTm';
ac_speech_training(p);

% score_speech(p); % Tom's script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end