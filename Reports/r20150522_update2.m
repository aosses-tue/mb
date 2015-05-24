function r20150522_update2(f)
% function r20150522_update2(f)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/05/2015
% Last update on: 22/05/2015 % Update this date manually
% Last use on   : 22/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

bDiary = 0;
Diary(mfilename,bDiary);
tic

% f1 = [1000 3000 5000];
% for i = 1:3
%     r20150522_update(f1(i));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params:
if nargin < 1
    f       = 1000;
end

stepJND     = [];
testJND     = [];
NJND        = [];
sigmaValues = [];
Nsigma      = [];
testLevels  = [];
Nlevels     = [];
Ntimes      = [];
Mat         = [];
CondN       = [];
CondCounter = [];

load([Get_TUe_paths('outputs') delim 'JNDcalculation-all-' num2str(f) '-'])

% JNDcalc    = nan(Nlevels,Nsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Creating test tone:
% fs      = 44100;
% dur     = 4; % in seconds
% y       = .5*Create_sin(f,dur,fs,0);
% rampup  = 5; % ms
% rampdn  = 5; % ms
% ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var2latex(Mat);

var2latex([sigmaValues' JNDcalc'])

toc

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

% Inline functions:
function outsig = Il_create_tone(insig,fs,lvl)

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(200e-3,fs); ...
               outsigtmp;
               Gen_silence(200e-3,fs)];
