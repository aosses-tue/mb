function r20150515_update
% function r20150515_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/05/2015
% Last update on: 18/05/2015 % Update this date manually
% Last use on   : 18/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc 

bDiary = 0;
Diary(mfilename,bDiary);
tic

stepJND     = .1;
testJND     = .7:stepJND:1.3;
% sigmaValues = [0 1 5 10]; % MU
sigmaValues = [.9 1 2]; % [0 0.1 0.5]; % MU
Nsigma      = length(sigmaValues);
% testLevels  = [20 40 60 80 100];
testLevels  = [60 80]; % [26 60 100];
Nlevels     = length(testLevels);
Ntimes      = 100; % 100 calculations for each condition

JNDcalc     = nan(Nlevels,Nsigma);
f           = 3000;
for i = 1:Nsigma
    
    opts = [];
    opts.f = f;
    opts.bPart3 = 1;
    opts.bPart4 = 0;
    opts.sigma  = sigmaValues(i);
    opts.sigmaTimes = Ntimes;
    opts.bDebug = 0;
    outs = r20150501_update_opt(opts);
    
    for j = 1:Nlevels
           
        opts = [];
        opts.bPart3     = 0;
        opts.bPart4     = 1;
        opts.f          = outs.f;
        opts.Criterion  = outs.Criterion;
        opts.testLevel  = testLevels(j);
        disp( sprintf('   - variance = %.2f: lvl = %.0f dB \n',sigmaValues(i),testLevels(j)) )
        opts.testJND    = testJND;
        opts.ytmp       = outs.ytmp;
        opts.crit       = outs.crit;
        opts.sigma      = sigmaValues(i);
        opts.sigmaTimes = Ntimes;
        opts.bDebug     = 1;
        outs2           = r20150501_update_opt(opts);
        
        JNDcalc(j,i)    = outs2.JNDcurrent;
        % JNDrecog(j,i)   = outs2.JNDrecognised;
    end
end

save([Get_TUe_paths('outputs') delim 'JNDcalculation'])
toc

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
