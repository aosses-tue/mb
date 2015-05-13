function y = r20150515_update(x)
% function y = r20150515_update(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/05/2015
% Last update on: 11/05/2015 % Update this date manually
% Last use on   : 11/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
tic

stepJND     = .05;
testJND     = .5:stepJND:1.8;
% sigmaValues = [1 5 10 20:2:30 40 50]; % MU
sigmaValues = [0 1 5 10]; % MU
Nsigma      = length(sigmaValues);
% testLevels  = [20 40 60 80 100];
testLevels  = [20 60 100];
Nlevels     = length(testLevels);
Ntimes      = 100; % 100 calculations for each condition

JNDcalc     = nan(Nlevels,Nsigma);

for i = 1:Nsigma
    
    opts = [];
    opts.bPart3 = 1;
    opts.bPart4 = 0;
    opts.sigma  = sigmaValues(i);
    opts.sigmaTimes = Ntimes;
    outs = r20150501_update_opt(opts);

    if i == 1
        disp( sprintf('variance = %.0f\n',opts.sigma) )
    end
    
    for j = 1:Nlevels
           
        opts = [];
        opts.bPart3     = 0;
        opts.bPart4     = 1;
        opts.f          = outs.f;
        opts.Criterion  = outs.Criterion;
        opts.testLevel  = testLevels(j);
        disp( sprintf('   - lvl = %.0f dB \n',testLevels(j)) )
        opts.testJND    = testJND;
        opts.ytmp       = outs.ytmp;
        opts.crit       = outs.crit;
        opts.sigma      = sigmaValues(i);
        opts.sigmaTimes = Ntimes;
        outs2           = r20150501_update_opt(opts);
        
        JNDcalc(j,i)    = outs2.JNDcurrent;
    end
end

save([Get_TUe_paths('outputs') delim 'JNDcalculation'])
toc

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
