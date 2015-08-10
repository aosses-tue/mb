function outs = r20150522_update_opt(opts)
% function outs = r20150522_update_opt(opts)
%
% 1. Description:
%       outs.JNDcurrent - corresponds to the estimated JND when a deviation
%                         of opts.sigma is considered.
% 
%       This script is being controlled from r20150522_update.m
% 
% 2. Stand-alone example:
%       r20150522_update; 
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20150522_update.m, r20150522_update2.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 18/05/2015
% Last update on: 22/05/2015 
% Last use on   : 22/05/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    opts = [];
end
outs = [];

opts = ef(opts,'bPart4',1);

opts = Ensure_field(opts,'sigma',.8);
opts = Ensure_field(opts,'sigmaTimes',1);
opts = Ensure_field(opts,'f',1000);
opts = Ensure_field(opts,'bDebug',0);

f = opts.f;

bPart4 = opts.bPart4;% Determining the noise deviation: Weber's approach
bDebug = opts.bDebug;

sigmaValues = opts.sigma;
sigmaTimes  = opts.sigmaTimes;

% Common parameters:
% close all
model   = 'dau1996';
fs      = 44100;
idx_compare = round(2*fs):round(3*fs);
Delta   = 1/fs;
t       = ( 1:44100 )/fs;
S       = ones(length(t),1);

if bPart4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same as 3, but now it looks for a threshold in levels other than 60 dB SPL:

if bDebug 
    Mat = [];
end

Criterion = opts.Criterion;

f           = opts.f;
testJND     = opts.testJND;

opts        = ef(opts,'testLevel',80);
testLevel   = opts.testLevel;
ytmp        = opts.ytmp;
crit        = opts.crit;
sigma       = sigmaValues;
testLevels  = testLevel*ones(size(testJND));

for idx = 1:length(testJND) 

    lvl1    = testLevels(idx);
    SMTc    = Il_adjust_tone(ytmp,fs,lvl1);
    
    lvl2    = subtract_dB( lvl1, testJND(idx) );
    SMTc2   = Il_adjust_tone(ytmp,fs,lvl2);
    
    if bDebug
        lvlTarget = sum_dB_arit([lvl1 lvl2]);
        Mat = [Mat; round(100*[lvlTarget lvl1 lvl2])/100];
        % var2latex(Mat);
    end
    
    absThres_at_f = absolutethreshold(f);
    
    if lvl2 < absThres_at_f
        warning('Test level below absolute threshold...')
        pause(5);
    end
    
    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    
    [RMTc   fc t] = Get_internal_representations_deterministic(SMTc ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);
    N = length(idx_compare);
    
    count_row = 1;

    for idx_times = 1:sigmaTimes
        
        mu  = 0;
        yn1 = normrnd(mu,sigma,N,1);
        yn2 = normrnd(mu,sigma,N,1);

        RMTc_n  = RMTc(idx_compare)  + yn1;
        RMTc2_n = RMTc2(idx_compare) + yn2;

        Mmin = floor(min(RMTc_n));
        Mmax = ceil(max(RMTc2_n));
        Centres = Mmin-1:1/16:Mmax+1;

        [PDF1 yi1 stats1] = Probability_density_function(RMTc_n,Centres);
        [PDF2 yi2 stats2] = Probability_density_function(RMTc2_n,Centres);
        
        [PDF1t yi1t stats1t] = Probability_density_function(RMTc(idx_compare),Centres);
        [PDF2t yi2t stats2t] = Probability_density_function(RMTc2(idx_compare),Centres);
        
        dprime_emp(count_row,idx) = ( stats2.mean - stats1.mean )/stats1.std;
        %dprime_emp2(count_row,idx) = ( stats2t.mean2 - stats1t.mean2 )/sigma;
        if idx_times == 1
            dprime(idx) = ( stats2t.mean - stats1t.mean )/sigma;
        end
        
        if bDebug
                        
            if idx_times == 1
                
                % figure; plot(yi1,PDF1), hold on; 
                % title(sprintf('Target level = %.1f dB, sigma=%.4f',sum_dB_arit([lvl1 lvl2]), sigma))
                % plot(yi2,PDF2,'r')
                % grid on
                % xlim([lvl1-15 lvl1+5]);
                % ylim([0 1])
                % 
                % fprintf('dprime = %.4f, noise normal = %s, noise+test normal = %s\n',dprime(count_row,idx),stats1.bNormal,stats2.bNormal);
                
            end
            fprintf('%.0f, noise1 = %.3f (bNormal %.0f), noise2 = %.3f (bNormal %.0f)\n', ...
                                                                        count_row, ...
                                                                        mean( yn1 ), ...
                                                                        Is_normal_distributed(yn1), ...
                                                                        mean( yn2 ), ...
                                                                        Is_normal_distributed(yn2));
        
        end
        count_row = count_row + 1;

    end
    
end

% dprime_not_met  = sum(dprime_emp<=Criterion);

try
    JND_3AFC = interp1(dprime,testJND,Criterion);
catch
    disp('interpolation failed, no JND found');
    JND_3AFC = nan(1);
end

outs.dprime         = dprime;
% outs.dprime_practice = dprime_emp;
% outs.JNDrecognised  = 100-dprime_not_met;
outs.JNDcurrent     = JND_3AFC;
outs.Mat            = Mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Inline functions:
function outsig = Il_adjust_tone(insig,fs,lvl)

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(200e-3,fs); ...
               outsigtmp;
               Gen_silence(200e-3,fs)];
           