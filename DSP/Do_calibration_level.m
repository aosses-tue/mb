function outsig = Do_calibration_level(method,insig,param)
% function outsig = Do_calibration_level(method,insig,param)
%
% 1. Description:
%       method / name   param
%       0 / AMT         desired level [dB SPL]
%       1 / Zwicker     audio file to calibrate (not required)
%       2 / AMT dB(A)   desired level [dB(A)]
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/08/2014
% Last update on: 21/08/2014 % Update this date manually
% Last use on   : 21/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

switch method
    case 0 % AMT calibration
        disp([mfilename '.m: AMT calibration (criterion: 0 dBFS = 100 dB SPL)'])
        fprintf('\t Your file will be adjusted to %.3f dBFS = %.3f dB SPL\n',param-100,param)
        outsig = setdbspl(insig,param);
    case 1 % Zwicker's calibration
        % NOT TESTED
        outsig =  insig;
    case 2 % AMT calibration considering A-weighting
        disp([mfilename '.m: AMT calibration (criterion: 0 dBFS = 100 dB(A))'])
        fprintf('\t Your file will be adjusted to %.3f dBFS = %.3f dB (A)\n',param-100,param)
        fs = 44100;
        warning('considering fs of 44100 Hz')
        [b a] = adsgn(44100);
        insigA = filter(b,a,insig);
        gain_in_dB = (param - 100) - rmsdb(insigA);
        gain = From_dB(gain_in_dB);
        % rmsdb(gain*insigA) % this is the value adjusted to lvl
        outsig = gain*insig;
        fprintf('\t Your file has been adjusted to an RMS of %.3f dBFS (Delta of %.3f dB)\n',rmsdb(outsig),gain)
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
