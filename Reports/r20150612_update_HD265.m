function r20150612_update_HD265
% function r20150612_update_HD265
%
% 1. Description:
%       Automatic generation of APEX experiments using the stimuli inside
%       the folder 'path'
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/06/2015
% Last update on: 11/06/2015 % Update this date manually
% Last use on   : 11/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bCreateAPEX = 1;

if bCreateAPEX
    path = 'D:\Documenten-TUe\02-Experiments\Set-up\calsignals-20150608\'; 
    savefilename = 'calibration-APEX-HD265-new.xml';
    opts.bCalTone = 1;
    opts.presentation = 'diotic';
    a3testacstimuli(savefilename, path,3,opts);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
