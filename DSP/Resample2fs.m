function y = Resample2fs(directory,fs)
% function y= Resample2fs(directory,fs)
%
% 1. Description:
%       To resample all wav files inside directory to a desired fs value
% 
% 2. Stand-alone example:
%       directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-calibrated\';
%       Resample2fs(directory);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 19/01/2015
% Last update on: 19/01/2015 
% Last use on   : 25/03/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

if nargin < 2
    fs = 44100;
end

if nargin == 0
    directory = uigetdir('Select a directory with wav-files');
end

file_orig = Get_filenames(directory,['*.wav']);
dirout = sprintf('%s%sresampled-at-%.0f-Hz%s',directory,delim,fs,delim);
Mkdir(dirout);

for i = 1:length(file_orig)
    filename(i) = file_orig(i);
    [x fs_old] = Wavread([directory delim filename{i}]);
    y = resample(x,fs,fs_old);
    Wavwrite(y,fs,[dirout filename{i}]);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
