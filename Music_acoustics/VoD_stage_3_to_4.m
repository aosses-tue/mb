function VoD_stage_3_to_4(directory)
% function VoD_stage_3_to_4(directory)
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

bDiary = 0;
Diary(mfilename,bDiary);

if nargin < 1
    directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test2\Stage3\';
end

acmode = [2 4];
field = 1:2;

for i = 1:length(acmode)
    mode = acmode(i)-1;
    for j = 1:length(field)
        basename = sprintf('%smode-%.0f-v_%.0f-',directory,mode,field(j));
        
        N = 4;
        [x1 fs] = Wavread([basename 's1r']);
        [x2 fs] = Wavread([basename 's1m']);
        [x3 fs] = Wavread([basename 's2r']);
        [x4 fs] = Wavread([basename 's2m']);
        
        yane = 1/N*( x1      + x3      );
        yrev = 1/N*( x1 + x2 + x3 + x4 );
        
        Wavwrite(yane,fs,[basename 'ane']);
        Wavwrite(yrev,fs,[basename 'rev']);
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
