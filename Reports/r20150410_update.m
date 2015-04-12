function y = r20150410_update(x)
% function y = r20150410_update(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 09/04/2015
% Last update on: 09/04/2015 % Update this date manually
% Last use on   : 09/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

path = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\JND-in-level\Stimuli\';

testref  = 'test-fc-5000-SPL-65dB.wav';
audioref = 'BBN-BW-9980-SPL-65dB.wav';

SPL     = 65;
[x fs] = Wavread([path audioref]);

SPLtest = [30 35 40];
att     = SPLtest-SPL;
bSaveWave = 1;

if bSaveWave
    
    fname = sprintf('%sBBN-BW-9980-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(1));
    Wavwrite(From_dB(att(1))*x,fs,fname)

    fname = sprintf('%sBBN-BW-9980-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(2));
    Wavwrite(From_dB(att(2))*x,fs,fname)
    
    fname = sprintf('%sBBN-BW-9980-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(3));
    Wavwrite(From_dB(att(3))*x,fs,fname)

end

testref  = 'test-fc-5000-SPL-65dB.wav';
SPL     = 65;
[x fs] = Wavread([path testref]);

SPLtest = [20 60 80];
att     = SPLtest-SPL;

if bSaveWave

    
    fname = sprintf('%stest-fc-5000-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(1));
    Wavwrite(From_dB(att(1))*x,fs,fname)

    fname = sprintf('%stest-fc-5000-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(2));
    Wavwrite(From_dB(att(2))*x,fs,fname)
    
    fname = sprintf('%stest-fc-5000-SPL-%.0fdB',Get_TUe_paths('outputs'),SPL+att(3));
    Wavwrite(From_dB(att(3))*x,fs,fname)

end



if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
