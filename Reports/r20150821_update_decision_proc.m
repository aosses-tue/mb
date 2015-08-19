function outs = r20150821_update_decision_proc(opts)
% function outs = r20150821_update_decision_proc(opts)
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
% Original file name: r20150522_update_opt.m
% Created on    : 18/08/2015
% Last update on: 18/08/2015 
% Last use on   : 18/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    opts = [];
end
outs = [];

opts = Ensure_field(opts,'sigma',.8);
opts = Ensure_field(opts,'sigmaTimes',1);
opts = Ensure_field(opts,'f',1000);
opts = Ensure_field(opts,'bDebug',0);
opts = Ensure_field(opts,'model','dau1996');
opts = Ensure_field(opts,'tmin',2);
opts = Ensure_field(opts,'tmax',3);
opts = Ensure_field(opts,'fs',44100);
opts = Ensure_field(opts,'bDeterministicTemplate',1);
opts = Ensure_field(opts,'siltime',200e-3);

bDeterministicTemplate = opts.bDeterministicTemplate;
siltime = opts.siltime;

f = opts.f;
tmin = opts.tmin;
tmax = opts.tmax;

bDebug = opts.bDebug;

sigmaValues = opts.sigma;
sigmaTimes  = opts.sigmaTimes;
model = opts.model;

% Common parameters:
fs      = opts.fs;
idx_compare = round(tmin*fs+1):round(tmax*fs);

if bDebug 
    Mat = [];
end

Criterion = opts.Criterion;

f           = opts.f;
testJND     = opts.testJND;

opts        = ef(opts,'testLevel',80);
testLevel   = opts.testLevel;
insig       = opts.insig;
crit        = opts.crit;
sigma       = sigmaValues;
testLevels  = testLevel*ones(size(testJND));

varout = [];
nCorrectAnswers = NaN;

for idx = 1:length(testJND) 

    if ( nCorrectAnswers(end) < 90 )| (isnan(nCorrectAnswers))
        
        if idx == 1 % lvl1 does not change
            lvl1= testLevels(idx);
            SMTc= Il_adjust_tone(insig,fs,lvl1,siltime);
        end

        lvl2    = subtract_dB( lvl1, testJND(idx) );
        SMTc2   = Il_adjust_tone(insig,fs,lvl2,siltime);

        if idx == 1 % then we generate the suprathreshold signal
            lvl3    = subtract_dB( lvl1, 5 );
            SMTc3   = Il_adjust_tone(insig,fs,lvl3,siltime);
        end

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

        if idx == 1
            [RM fc t] = Get_internal_representations_deterministic(SMTc,fs,model,setup);
        end
        
        [RMTc  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);

        if idx == 1
            
            RMTsupra = Get_internal_representations_deterministic(SMTc+SMTc3,fs,model,setup); % rmsdb(SMTc+SMTc3)+100
            
            [N M] = size(RMTsupra);
            
            if bDeterministicTemplate == 0
                
                mu = 0;
                RM_avg       = il_get_avg_ir(RM(idx_compare,:)      ,mu,sigma,sigmaTimes);
                RMTsupra_avg = il_get_avg_ir(RMTsupra(idx_compare,:),mu,sigma,sigmaTimes);
                                
            else
                RM_avg       = RM;
                RMTsupra_avg = RMTsupra;
            end
            
            T = Get_template_append(RM_avg,RMTsupra_avg,fs);
            N = length(idx_compare);
            M = size(T,2);

        end

        for k = 1:sigmaTimes
            mu  = 0;
            yn1 = normrnd(mu,sigma,N,M);
            yn2 = normrnd(mu,sigma,N,M);

            yn3 = normrnd(mu,sigma,N,M);
            yn4 = normrnd(mu,sigma,N,M);
            yn5 = normrnd(mu,sigma,N,M);

            RMTc_n  = RM(idx_compare,:)  + yn1;
            RMTc2_n = RMTc(idx_compare,:) + yn2;

            RMTc_n3 = RM(idx_compare,:) + yn3;
            RMTc_n4 = RM(idx_compare,:) + yn4;
            RMTc_n5 = RM(idx_compare,:) + yn5;

            Mmin = floor(min(RMTc_n));
            Mmax = ceil(max(RMTc2_n));
            
            tmp_mue1 = optimaldetector(RMTc_n - RMTc_n3,T); % Noise alone
            tmp_mue2 = optimaldetector(RMTc_n - RMTc_n4,T); % Noise alone
            [mue(1,idx) mue(2,idx)] = max([tmp_mue1 tmp_mue2]);
            mue(3,idx) = optimaldetector(RMTc2_n - RMTc_n5,T); % Signal + noise
            mue(4,idx) = mue(3,idx) - mue(1,idx);

            tmp_mue1 = optimaldetector( RMTc_n,T);
            tmp_mue2 = optimaldetector(RMTc_n4,T);
            [mue2(1,idx) mue2(2,idx)] = max([tmp_mue1 tmp_mue2]);
            mue2(3,idx) = optimaldetector(RMTc2_n,T);
            mue2(4,idx) = mue(3,idx) - mue(1,idx);

            dprime(idx)     = mue(4,idx)/sigma;
            varout(k,idx)   = mue(4,idx)/sigma;
            varout2(k,idx)  = mue2(4,idx)/sigma;
            nCorrectAnswers = sum(varout > opts.Criterion);
        end
    else
        
        dprime(idx) = NaN;
        varout(k,idx) = NaN;
        varout2(k,idx) = NaN;
        nCorrectAnswers = NaN;
        
    end
    
end

try
    idx = find(dprime > 1); % to use less numbers in the interpolation
    JND_3AFC = interp1(dprime(idx),testJND(idx),Criterion);
catch
    disp('interpolation failed, no JND found');
    JND_3AFC = nan(1);
end

outs.dprime         = dprime;
outs.JNDcurrent     = JND_3AFC;
outs.Mat            = Mat;
outs.varout = varout;
outs.varout2 = varout2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inline functions:
function outsig = Il_adjust_tone(insig,fs,lvl,siltime)

if nargin <4
    siltime = 200e-3;
end

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(siltime,fs); ...
               outsigtmp;
               Gen_silence(siltime,fs)];
           
function RMTout = il_get_avg_ir(RMTin,mu,sigma,sigmaTimes)

[N M] = size(RMTin);
Rtmp = [];

for k = 1:sigmaTimes
    yn1 = normrnd(mu,sigma,N,M);
    Rtmp = [Rtmp RMTin(:)+yn1(:)];
end
Rtmp = mean(Rtmp,2);
RMTout = reshape(Rtmp,N,M);