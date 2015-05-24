function outs = r20150501_update_opt(opts)
% function outs = r20150501_update_opt(opts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/05/2015
% Last update on: 12/05/2015 % Update this date manually
% Last use on   : 12/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    opts = [];
end
outs = [];

opts = ef(opts,'bPart3',1);
opts = ef(opts,'bPart4',1);

opts = Ensure_field(opts,'sigma',30);
opts = Ensure_field(opts,'sigmaTimes',100);
opts = Ensure_field(opts,'f',1000);
opts = Ensure_field(opts,'bDebug',0);

f = opts.f;

bPart3 = opts.bPart3; % Determining the noise deviation: Weber's approach
bPart4 = opts.bPart4;
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

if bPart3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes respect bPart2:
% 
%   - Estimate threshold for N times
%   - No suprathreshold signal
%   - No template used

testJND     = 1.0; 
testLevels  = 60*ones(size(testJND)); % Level test tone

dur     = 4; % in seconds

y       = .5*Create_sin(f,dur,fs,0);
rampup  = 5; % ms
rampdn  = 5; % ms
ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);

durtot  = dur+400e-3;

for idx = 1:length(testJND) 

    lvl1    = testLevels(idx);
    SMTc    = Il_create_tone(ytmp,fs,lvl1);
    
    lvl2    = subtract_dB( lvl1, testJND(idx) );
    SMTc2   = Il_create_tone(ytmp,fs,lvl2);
    
    if bDebug
        lvlTarget = sum_dB_arit([lvl1 lvl2]);
        Mat = round(100*[lvlTarget lvl1 lvl2])/100;
        var2latex(Mat);
    end
    
    % sound(SMTc ,fs)
    % sound(SMTc2,fs)

    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    
    [RMTc   fc t] = Get_internal_representations_deterministic(SMTc      ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);
    N = length(RMTc);
    
    count_row = 1;

    mu = 0;
    sigma = sigmaValues;

    for idx_times = 1:sigmaTimes
        yn1 = normrnd(mu,sigma,N,1); 
        yn2 = normrnd(mu,sigma,N,1);

        RMTc_n  = RMTc  + yn1;
        RMTc2_n = RMTc2 + yn2;

        diffVector                  = RMTc2_n(idx_compare) - RMTc_n(idx_compare);
        difference(count_row,idx) = 1/fs*sum(diffVector); % current difference
        
        if bDebug
            figure; plot(diffVector);
            
            %M = 1000;
            Mmin = floor(min(RMTc_n));
            Mmax = ceil(max(RMTc2_n));
            
            Centres = Mmin-1:Mmax+1;
            [PDF1 yi1] = Probability_density_function(RMTc_n,Centres);
            [PDF2 yi2] = Probability_density_function(RMTc2_n,Centres);
            
            figure; 
            plot(yi1,PDF1, yi2,PDF2); grid on 
            xlim([lvl1-15 lvl1+15]);
            ylim([0 1]);
        end
        count_row = count_row + 1;

    end
    
end

M0dB = mean( difference(:,end) );
% difference = difference/M0dB;

idCum = cumsum(sigmaTimes);

crit    = 70.7; % percentage correct
P3      = percentile(difference,crit);

Criterion = P3; % interp1(testJND,P3,0.95);

outs.f  = f;
outs.PC = P3;
outs.Criterion = Criterion;
outs.crit   = crit;
outs.ytmp = ytmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same as 3, but now it looks for a threshold in levels other than 60 dB SPL:

if ~bPart3
    if ~isfield(opts,'Criterion')
        error('Set bPart3 first');
    else
        Criterion = opts.Criterion;
    end
end

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
    SMTc    = Il_create_tone(ytmp,fs,lvl1);
    
    lvl2    = subtract_dB( lvl1, testJND(idx) );
    SMTc2   = Il_create_tone(ytmp,fs,lvl2);
    
    if bDebug
        lvlTarget = sum_dB_arit([lvl1 lvl2]);
        Mat = round(100*[lvlTarget lvl1 lvl2])/100;
        var2latex(Mat);
    end
    
    absThres_at_f = absolutethreshold(f);
    
    if lvl2 < absThres_at_f
        warning('Test level below absolute threshold...')
        pause(5);
    end
    
    % sound(SMTc ,fs)
    % sound(SMTc2,fs)

    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    
    [RMTc   fc t] = Get_internal_representations_deterministic(SMTc ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);
    N = length(RMTc);
    
    count_row = 1;

    for idx_times = 1:sigmaTimes
        
        mu = 0;
        yn1 = normrnd(mu,sigma,N,1);
        yn2 = normrnd(mu,sigma,N,1);

        RMTc_n  = RMTc  + yn1;
        RMTc2_n = RMTc2 + yn2;

        diffVector                = RMTc2_n(idx_compare) - RMTc_n(idx_compare);
        diffTest(count_row,idx)   = 1/fs*sum(diffVector); % current difference
        
        if bDebug
                        
            if idx_times == 1
                
                % figure; plot(diffVector);
                
                %M = 1000;
                Mmin = floor(min(RMTc_n));
                Mmax = ceil(max(RMTc2_n));

                Centres = Mmin-1:1/16:Mmax+1;

                [PDF1 yi1 stats1] = Probability_density_function(RMTc_n(idx_compare),Centres);
                figure; plot(yi1,PDF1), hold on; 
                title(sprintf('Target level = %.1f dB, sigma=%.4f',sum_dB_arit([lvl1 lvl2]), sigma))
                [PDF2 yi2 stats2] = Probability_density_function(RMTc2_n(idx_compare),Centres);
                plot(yi2,PDF2,'r')
                grid on
                xlim([lvl1-15 lvl1+5]);
                ylim([0 1])
                dprime = ( stats2.mean - stats1.mean )/stats1.std;
                
                fprintf('dprime = %.4f, noise normal = %s, noise+test normal = %s\n',dprime,stats1.bNormal,stats2.bNormal);
            end
        
        end
        count_row = count_row + 1;

    end
    
end

M0dB = mean( diffTest(:,end) );
P3new  = percentile(diffTest,crit);

amount = sum(diffTest<=Criterion);

% to avoid problems with interpolation...
% idx2delete = find(diff(P3new)==0);
% P3new(idx2delete) = [];
% testJND(idx2delete) = [];
 
tmp = interp1(P3new,testJND,Criterion);
% tmp = P3new;
outs.JNDrecognised = 100-amount;
outs.JNDcurrent = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Inline functions:
function outsig = Il_create_tone(insig,fs,lvl)

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(200e-3,fs); ...
               outsigtmp;
               Gen_silence(200e-3,fs)];
           