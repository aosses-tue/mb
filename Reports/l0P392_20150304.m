function l0P392_20150304
% function l0P392_20150304
%
% 1. Description:
%	For assignment of the advanced perception course (after Speech 
%   perception 2)
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 02/03/2015
% Last update on: 02/03/2015 % Update this date manually
% Last use on   : 02/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bFigures4lesson = 1;

dir = 'D:\MATLAB_git\Reports\l0P392_20150304\';

% dirNL = [dir 'Dutch.wav\'];
% dirEN = [dir 'English.wav\'];
dirEN = [dir 'test\'];

bEngels = 1;

cd(dir);

% files = Get_filenames(dirNL,'*.wav');
if bEngels
    dir2use = dirEN;
else
    dir2use = dirNL;
end

files = Get_filenames(dir2use,'*.wav');
    
B = 3:10;
rate = 16:-1:1;

opts.bSave = 1;
opts.bPlot = 0;

bSmearBits = 1;
bSmearFs = 1;

for i = 1:length(files)
    
    if bSmearBits
        for j = 1:length(B)
            d_quant([dir2use files{i}],B(j),opts);
        end
    end
    
    if bSmearFs
        for j = 1:length(rate)
            d_resample([dir2use files{i}],rate(j));
        end
    end
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
