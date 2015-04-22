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

bPart1 = 0;
bPart2 = 1;

if bPart1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart2
    
    fs = 44100;
    N = 200e-3*44100;
    sig = [1; zeros(N-1,1)];
    f_abt = 1/2e-3; % 2-ms sampling period
    
    % Calculation of coefficients of critical band filterbank
    S = ch_make_fttbank1(fs);
    
    % Applying critical band filterbank
    [fgrp, S] = ch_ftt_bank1(sig, S, f_abt,fs);
    
    for chan = 1:24
        CBvalue(:,chan) = S.y4{chan};
    end
    
    figure;
    freqz(CBvalue(:,5),1,4096)
    disp('')
    
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
