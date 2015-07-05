function y = r20150608_test_stimuli
% function y = r20150608_test_stimuli
%
% 1. Description:
%       It generates the test stimuli to be used in the periodic verification
%       (calibration) of the auditoryLab
%
% 2. Stand-alone example:
%       r20150608_test_stimuli;
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 08/06/2015
% Last update on: 08/06/2015 % Update this date manually
% Last use on   : 08/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
outdir = 'D:\Documenten-TUe\02-Experiments\Set-up\calsignals\';

f = [125 250 500 1000 2000 4000 8000];
fs = 44100;

for i = 1:length(f)
    y = Create_sin(f(i),10,fs);
    Wavwrite(    y,fs,[outdir 'tone-' num2str(f(i)) 'Hz-Amp-100.wav']); %   0 dB
    Wavwrite(0.5*y,fs,[outdir 'tone-' num2str(f(i)) 'Hz-Amp-050.wav']); %  -6 dB
    Wavwrite(0.1*y,fs,[outdir 'tone-' num2str(f(i)) 'Hz-Amp-010.wav']); % -20 dB
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
