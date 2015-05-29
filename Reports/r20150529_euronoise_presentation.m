function r20150529_euronoise_presentation
% function r20150529_euronoise_presentation
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/05/2015
% Last update on: 28/05/2015 % Update this date manually
% Last use on   : 28/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bRealignSignals = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bRealignSignals

dirp = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\2015-02-wav-files\02-calibrated\';
f1 = [dirp 'model-ac-2-dist-ane'];
f2 = [dirp 'model-ac-2-dist-rev'];
f3 = [dirp 'model-ac-4-dist-ane'];
f4 = [dirp 'model-ac-4-dist-rev'];

fs = 44100;
to_crop_predict = round( [0 0.06 0.045 0]*fs )+1;

for i = 4 %1:4
    N       = to_crop_predict(i);
    exp1    = sprintf('[x fs] = Wavread([f%.0f ''.wav'']);',i);
    eval(exp1);
    y = x(N:end);
    
    exp1 = sprintf('Wavwrite(y,fs,[f%.0f ''-aligned.wav'']);',i);
    eval(exp1);
end

dirm = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\2015-02-wav-files\02-calibrated\';
f1 = [dirm 'meas-ac-2-dist-ane'];
f2 = [dirm 'meas-ac-2-dist-rev'];
f3 = [dirm 'meas-ac-4-dist-ane-HP'];
f4 = [dirm 'meas-ac-4-dist-rev-HP'];

to_crop_meas = round( [0.23 0 0 0.04]*fs )+1; 

for i = 4 % 1:4
    N       = to_crop_meas(i);
    exp1    = sprintf('[x fs] = Wavread([f%.0f ''.wav'']);',i);
    eval(exp1);
    y = x(N:end);
    
    exp1 = sprintf('Wavwrite(y,fs,[f%.0f ''-aligned.wav'']);',i);
    eval(exp1);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
