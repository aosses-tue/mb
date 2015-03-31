function [outs] = demo_adaptloop2_gen_stim(opts)
% function [outs] = demo_adaptloop2_gen_stim(opts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 24/03/2015
% Last update on: 24/03/2015 % Update this date manually
% Last use on   : 24/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

if nargin < 1
    opts = [];
end

opts = Ensure_field(opts,'outputs',Get_TUe_paths('outputs'));
opts = Ensure_field(opts,'f',1000);

f   = opts.f;
dur = 200e-3;
dur_sil = 50e-3;

fs  = 44100;
y   = Create_sin(f,dur,fs,0);
SPL = 65;
y   = setdbspl(y,SPL);
y   = [y; Gen_silence(dur_sil,fs)];
fname = sprintf('%.0f-sine-%.0f-dBSPL',f,SPL);
fname = [opts.outputs fname];

Wavwrite(y,fs,fname);

outs.fname = [fname '.wav'];

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
