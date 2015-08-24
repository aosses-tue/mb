function outs = r20150821_update_decision_proc_exp(insigN,insigS,opts)
% function outs = r20150821_update_decision_proc_exp(insigN,insigS,opts)
%
% 1. Description:
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20150522_update.m, r20150522_update2.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 21/08/2015
% Last update on: 21/08/2015 
% Last use on   : 21/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    opts = [];
end
outs = [];

opts = Ensure_field(opts,'sigmaTimes',1);
opts = Ensure_field(opts,'f',1000);
opts = Ensure_field(opts,'bDebug',0);
opts = Ensure_field(opts,'model','dau1996');
opts = Ensure_field(opts,'tmin',2);
opts = Ensure_field(opts,'tmax',3);
opts = Ensure_field(opts,'fs',44100);
opts = Ensure_field(opts,'bDeterministicTemplate',1);
opts = Ensure_field(opts,'siltime',200e-3);

model = opts.model;

switch model
    case 'dau1996'
        opts = Ensure_field(opts,'sigma',.8);
    case 'jepsen2008'
        opts = Ensure_field(opts,'sigma',.85);
end

bDeterministicTemplate = 0;
siltime = opts.siltime; % maskers

f = opts.f;

opts = Ensure_field(opts,'fmin',f);
opts = Ensure_field(opts,'fmax',f);

fmin = opts.fmin;
fmax = opts.fmax;
tmin = opts.tmin;
tmax = opts.tmax;
ramp = opts.ramp;
gainN = opts.gainN;
AboveThres = opts.AboveThres;
bDebug = opts.bDebug;

sigma = opts.sigma;
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
testlevels  = opts.testlevels;

opts        = ef(opts,'testLevel',80);
testLevel   = opts.testLevel;
% insig       = opts.insig;
% crit        = opts.crit;
testLevels  = testLevel*ones(size(testlevels));
bInsig_is_buffer = 1;

varout = [];
nCorrectAnswers = NaN;

dprime = 100;
insigN_orig = insigN;
insigS_orig = insigS;
idx = 1;
currentlvl = AboveThres;
cont = 1;
nReversals = 0;

Staircase = [];
Reversals = [];
Step = 12;
bCorrect = 1;

while nReversals <= 5

    if bInsig_is_buffer

        if idx == 1
            insigsupra = insigS(1:tmax*fs);
        end
        
        insigtmp = il_randomise_insig(insigN_orig,fs,ramp);
        insigN = insigtmp(1:tmax*fs);
        insigtmp = il_randomise_insig(insigN_orig,fs,ramp);
        insigN2 = insigtmp(1:tmax*fs);
        insigtmp = il_randomise_insig(insigN_orig,fs,ramp);
        insigN3 = insigtmp(1:tmax*fs);
        insigtmp = il_randomise_insig(insigN_orig,fs,ramp);
        insigN4 = insigtmp(1:tmax*fs);
        insigtmp = il_randomise_insig(insigN_orig,fs,ramp);
        insigN5 = insigtmp(1:tmax*fs);
    end
        
    SMTN1   = gaindb(insigN ,gainN); % Masker
    SMTN2   = gaindb(insigN2,gainN); % Masker
    SMTN3   = gaindb(insigN3,gainN); % Masker
    SMTN4   = gaindb(insigN4,gainN); % Masker
    SMTN5   = gaindb(insigN5,gainN); % Masker
    SMTc2   = gaindb(insigS,currentlvl);
    SMTc3   = gaindb(insigsupra,AboveThres);
        
    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    setup.fmin      = fmin;
    setup.fmax      = fmax;

    [RM  fc t] = Get_internal_representations_deterministic(SMTN1,fs,model,setup);
    [RM2 fc t] = Get_internal_representations_deterministic(SMTN2,fs,model,setup);
    [RM3 fc t] = Get_internal_representations_deterministic(SMTN3,fs,model,setup);
    [RM4 fc t] = Get_internal_representations_deterministic(SMTN4,fs,model,setup);
    [RM5 fc t] = Get_internal_representations_deterministic(SMTN5,fs,model,setup);
    
    if size(RM,1) ~= size(SMTN1,1);
        if idx == 1
            warning('Downsampled internal representation')
            factor = round(size(SMTN1,1)/size(RM,1) );
            idx_compare = round(idx_compare / factor);
            if idx_compare(1) == 0
                idx_compare(1) = [];
            end
            idx_compare = idx_compare(1:factor:end);
            idxtodelete = length(idx_compare) - size(RM,1);
            idx_compare(end-idxtodelete+1) = [];
        end
    end
        
    Input2 = SMTN1+SMTc2; % sound(Input2,fs);
    [RMTc  fc t] = Get_internal_representations_deterministic(Input2,fs,model,setup);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To get the template:
    if idx == 1
        % InputSupra = SMTN1+SMTc3; % sound(InputSupra)
        
        insigN_tmp = gaindb(insigN_orig,gainN);
        [RM_avg, RMTsupra_avg, fs_intrep, outtmp] = Get_internalrep_stochastic(insigN_tmp, SMTc3,fs,model,sigma); %,Ntimes,idx_fc);
        if strcmp(model,'dau1996')
            idx_tmp = find(outtmp.fc >= fc(1) & outtmp.fc <= fc(end));
            RM_avg = RM_avg(idx_compare,idx_tmp);
            RMTsupra_avg = RMTsupra_avg(idx_compare,idx_tmp);
        end
        
        % RMTsupra = Get_internal_representations_deterministic(InputSupra,fs,model,setup); % rmsdb(SMTc+SMTc3)+100
        idx = 0;
        % [N M] = size(RMTsupra);
        % 
        % if bDeterministicTemplate == 0
        % 
        %     mu = 0;
        %     RM_avg       = il_get_avg_ir(RM(idx_compare,:)      ,mu,sigma,sigmaTimes);
        %     RMTsupra_avg = il_get_avg_ir(RMTsupra(idx_compare,:),mu,sigma,sigmaTimes);
        % 
        % else
        %     RM_avg       = RM;
        % end

        T = Get_template_append(RM_avg,RMTsupra_avg,fs);
        N = length(idx_compare);
        M = size(T,2);

    end
    % end find template
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    dprimetmp = [];
    for k = 1:sigmaTimes
        mu  = 0;
        yn1 = normrnd(mu,sigma,N,M);
        RM1_n = RM(idx_compare,:) + yn1;
        
        yn1 = normrnd(mu,sigma,N,M);
        RM2_n  = RM2(idx_compare,:) + yn1;
        
        yn1 = normrnd(mu,sigma,N,M);
        RM3_n  = RM3(idx_compare,:) + yn1;
        
        yn1 = normrnd(mu,sigma,N,M);
        RM4_n  = RM4(idx_compare,:) + yn1;
        
        yn1 = normrnd(mu,sigma,N,M);
        RM5_n  = RM5(idx_compare,:) + yn1;
        
        yn1 = normrnd(mu,sigma,N,M);
        RMTc_n = RMTc(idx_compare,:) + yn1;

        idx_com = idx_compare; % round(2*fs+1):round(3*fs);
        tmp_mue1 = optimaldetector(RM1_n(idx_com,:) - RM2_n(idx_com,:),T(idx_com,:)); % Noise alone
        tmp_mue2 = optimaldetector(RM3_n(idx_com,:) - RM4_n(idx_com,:),T(idx_com,:)); % Noise alone
        [mue(1,cont) mue(2,cont)] = max([tmp_mue1 tmp_mue2]); % noise sample with larger dprime is chosen
        
        mue(3,cont) = optimaldetector(RMTc_n - RM5_n,T); % Signal + noise
        mue(4,cont) = mue(3,cont) - mue(1,cont);

        dprimetmp(end+1)     = mue(4,cont)/sigma;
        cont = cont + 1;
    end
    
    dprime(end+1) = mean(dprimetmp);
        
    Staircase = [Staircase currentlvl];
    if dprime(end) > Criterion
        currentlvl = currentlvl - Step;
        bCorrect = 1;
        bWrong = 0;
    else
        if bCorrect == 1
            % bWrong = 1;
            bCorrect = 0;
        else
            currentlvl = currentlvl + Step;
            nReversals = nReversals + 1;

            if mod(nReversals,2) == 0

                Step = Step/2;

            end

            Reversals = [Reversals currentlvl];
            bCorrect = 1;
        end
    end
        
    if dprime(end) < 1.26
        disp('')
    else
        disp(num2str(cont))
    end
end

JND_3AFC = median(Reversals(end-round(nReversals/2):end));

outs.dprime         = dprime(end);
outs.JNDcurrent     = JND_3AFC;
outs.Reversals      = Reversals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inline functions:
function outsig = Il_adjust_tone(insig,fs,lvl,siltime,ramp)

if nargin <5
    ramp = 0;
end

if nargin <4
    siltime = 200e-3;
end

outsigtmp = Do_cos_ramp(insig,fs,ramp);

outsigtmp   = setdbspl(outsigtmp,lvl);

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

function outsig = il_randomise_insig(insig,fs,ramp)

if nargin < 3
    ramp = 0;
end

Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];

outsig = Do_cos_ramp(outsig,fs,ramp);
