function VoD_stage_3_to_4(directory,GainS1,GainS2)
% function VoD_stage_3_to_4(directory,GainS1,GainS2)
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
% Last use on   : 23/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

if nargin < 1
    % directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-test2\Stage3\';
    % directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm\Stage3\'; % including Harmonic
    dirmain   = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm-G1\'; % including Harmonic
end
dirin       = [dirmain 'Stage3' delim];
dirout      = [dirmain 'Stage4' delim];

if nargin < 2
    GainS1 = 2;
end

if nargin < 3
    GainS2 = 0;
end

sprintf('Gain to be applied to inlet signals S1: %.2f [dB]\n',GainS1)
sprintf('Gain to be applied to outlet signals S2: %.2f [dB]\n',GainS2)

acmode = [2 4];
field = 1:2;

for i = 1:length(acmode)
    mode = acmode(i)-1;
    for j = 1:length(field)
        basename = sprintf('%smode-%.0f-v_%.0f-',dirin,mode,field(j));
        baseout  = sprintf('%smode-%.0f-v_%.0f-',dirout,mode,field(j));
        
        N = 4;
        [x1 fs] = Wavread([basename 's1r']);
        [x2 fs] = Wavread([basename 's1m']);
        [x3 fs] = Wavread([basename 's2r']);
        [x4 fs] = Wavread([basename 's2m']);
        
        yane = 1/N*( From_dB(GainS1)* x1       + From_dB(GainS2)* x3      );
        yrev = 1/N*( From_dB(GainS1)*(x1 + x2) + From_dB(GainS2)*(x3 + x4));
        
        Wavwrite(yane,fs,[baseout 'ane']);
        Wavwrite(yrev,fs,[baseout 'rev']);
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
