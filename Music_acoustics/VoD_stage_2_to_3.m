function VoD_stage_2_to_3(directory)
% function VoD_stage_2_to_3(directory)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/02/2015
% Last update on: 13/02/2015 % Update this date manually
% Last use on   : 13/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

if nargin < 1
    % directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test-all-sources\Stage2\';
    directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm\Stage2\';
end

opts.bExtension = 1;
filenames = Get_filenames(directory,'*.txt',opts);

for i = 1:length(filenames)
    files = [directory filenames{i}];
    VoD_one_predicted_file_from_txt(files);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
