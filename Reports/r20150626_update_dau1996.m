function r20150626_update_dau1996(f, testLevels)
% function r20150626_update_dau1996(f, testLevels)
%
% 1. Description:
%
% 2. Stand-alone example:
%       f = 1000;
%       testLevels = 60;
%       r20150626_update_dau1996(f, testLevels);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 26/06/2015
% Last update on: 26/06/2015 % Update this date manually
% Last use on   : 26/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

bDiary = 0;
Diary(mfilename,bDiary);
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params:
if nargin < 1
    f       = 5000;
end

if nargin < 2
    testLevels  = [20 26 40 60 80 100];
end

testJND     = 0.8;
NJND        = length(testJND);
sigmaValues = [.3  .5:.1:1]; % [0 0.1 0.5]; % MU
Nsigma      = length(sigmaValues);

Nlevels     = length(testLevels);
Ntimes      = 1; % 1 calculations for each condition
JNDcalc     = nan(Nlevels,Nsigma);
Mat         = nan(NJND*Nlevels,3+Nsigma);
CondN       = Nsigma * Nlevels;
CondCounter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating test tone:
fs      = 44100;
dur     = 4; % in seconds
y       = .5*Create_sin(f,dur,fs,0);
rampup  = 5; % ms
rampdn  = 5; % ms
ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Nsigma

    crit            = 70.7;
    Criterion       = 1.26; % D-prime for 3-AFC    
    
    for j = 1:Nlevels % each testLevel
           
        opts = [];
        opts.bPart3     = 0;
        opts.bPart4     = 1;
        opts.f          = f;
        opts.Criterion  = Criterion;
        opts.testLevel  = testLevels(j);
        disp( sprintf('   - variance = %.2f: lvl = %.0f dB \n',sigmaValues(i),testLevels(j)) );
        opts.testJND    = testJND;
        opts.ytmp       = ytmp;
        opts.crit       = crit;
        opts.sigma      = sigmaValues(i);
        opts.sigmaTimes = Ntimes;
        opts.bDebug     = 1;
        outs2           = r20150522_update_opt(opts);
        
        Mat(1 + (j-1)*NJND:j*NJND,1:3) = outs2.Mat;
        Mat(1 + (j-1)*NJND:j*NJND,3+i) = transpose(outs2.dprime);
        
        JNDcalc(j,i)    = outs2.JNDcurrent;
        fprintf('Completed %.2f%% (%.0f out of %.0f conditions)\n',CondCounter/CondN*100,CondCounter,CondN);
        CondCounter     = CondCounter + 1;
        % JNDrecog(j,i)   = outs2.JNDrecognised;
    end
end

var2latex(Mat);

toc
save([Get_TUe_paths('outputs') delim 'JNDcalculation-all-' num2str(f) '-'])

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
