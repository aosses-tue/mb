function out = Adjust_speech_rms_proc(p, wav)
% function out = Adjust_speech_rms_proc(p, wav)
%
% Adjusts rms: some kind of calibration
%
% Programmed by Matthias Milczynski, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case 0
        out = feval(mfilename, []);
    case 1
        p = Ensure_field(p, 'rms_gain_dB', 0);
        out = p;
    case 2
        out = eq_rms(wav, rmsdb(wav) + p.rms_gain_dB);
        
end
