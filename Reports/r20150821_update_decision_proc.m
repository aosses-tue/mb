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

f = opts.f;
tmin = opts.tmin;
tmax = opts.tmax;

bDebug = opts.bDebug;

sigmaValues = opts.sigma;
sigmaTimes  = opts.sigmaTimes;
model = opts.model;

% Common parameters:
fs      = 44100;
idx_compare = round(tmin*fs+1):round(tmax*fs);
Delta   = 1/fs;
t       = ( 1:44100 )/fs;
S       = ones(length(t),1);

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

varout = [];
nCorrectAnswers = NaN;

for idx = 1:length(testJND) 

    if ( nCorrectAnswers(end) < 90 )| (isnan(nCorrectAnswers))
        lvl1    = testLevels(idx);
        SMTc    = Il_adjust_tone(ytmp,fs,lvl1);

        lvl2    = subtract_dB( lvl1, testJND(idx) );
        SMTc2   = Il_adjust_tone(ytmp,fs,lvl2);

        if idx == 1
            lvl3    = subtract_dB( lvl1, 5 );
            SMTc3   = Il_adjust_tone(ytmp,fs,lvl3);
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

        [RMTc   fc t] = Get_internal_representations_deterministic(SMTc ,fs,model,setup);
        [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);

        if idx == 1
            RMTsupra = Get_internal_representations_deterministic(SMTc+SMTc3,fs,model,setup);
            T = Get_template(RMTc,RMTsupra,fs);

            [xx xx Texcerpt] = Normalise_signal(T(idx_compare),fs);

        end

        N = length(idx_compare);

        count_row = 1;

        for k = 1:sigmaTimes
            mu  = 0;
            yn1 = normrnd(mu,sigma,N,1);
            yn2 = normrnd(mu,sigma,N,1);

            yn3 = normrnd(mu,sigma,N,1);
            yn4 = normrnd(mu,sigma,N,1);
            yn5 = normrnd(mu,sigma,N,1);

            RMTc_n  = RMTc(idx_compare)  + yn1;
            RMTc2_n = RMTc2(idx_compare) + yn2;

            RMTc_n3 = RMTc(idx_compare) + yn3;
            RMTc_n4 = RMTc(idx_compare) + yn4;
            RMTc_n5 = RMTc(idx_compare) + yn5;

            Mmin = floor(min(RMTc_n));
            Mmax = ceil(max(RMTc2_n));
            Centres = Mmin-1:1/16:Mmax+1;

            [PDF1 yi1 stats1] = Probability_density_function(RMTc_n,Centres);
            [PDF2 yi2 stats2] = Probability_density_function(RMTc2_n,Centres);

            [PDF1t yi1t stats1t] = Probability_density_function(RMTc(idx_compare),Centres);
            [PDF2t yi2t stats2t] = Probability_density_function(RMTc2(idx_compare),Centres);

            % mue(1,idx) = optimaldetector(RMTc_n - RMTc_n3,T(idx_compare));
            % mue(2,idx) = optimaldetector(RMTc_n - RMTc_n4,T(idx_compare));
            % mue(3,idx) = optimaldetector(RMTc2_n - RMTc_n,T(idx_compare));
            % 
            % mue2(1,idx) = optimaldetector( RMTc_n,T(idx_compare));
            % mue2(2,idx) = optimaldetector(RMTc_n3,T(idx_compare));
            % mue2(3,idx) = optimaldetector(RMTc2_n,T(idx_compare));

            tmp_mue1 = optimaldetector(RMTc_n - RMTc_n3,Texcerpt); % Noise alone
            tmp_mue2 = optimaldetector(RMTc_n - RMTc_n4,Texcerpt); % Noise alone
            [mue(1,idx) mue(2,idx)] = max([tmp_mue1 tmp_mue2]);
            mue(3,idx) = optimaldetector(RMTc2_n - RMTc_n5,Texcerpt); % Signal + noise
            mue(4,idx) = mue(3,idx) - mue(1,idx);

            tmp_mue1 = optimaldetector( RMTc_n,Texcerpt);
            tmp_mue2 = optimaldetector(RMTc_n4,Texcerpt);
            [mue2(1,idx) mue2(2,idx)] = max([tmp_mue1 tmp_mue2]);
            mue2(3,idx) = optimaldetector(RMTc2_n,Texcerpt);
            mue2(4,idx) = mue(3,idx) - mue(1,idx);

            dprime(idx) = ( stats2t.mean - stats1t.mean )/sigma;

            varout(k,idx) = mue(4,idx)/sigma;
            varout2(k,idx) = mue2(4,idx)/sigma;
            nCorrectAnswers = sum(varout > opts.Criterion);
        end
    else
        
        dprime(idx) = NaN;
        varout(k,idx) = NaN;
        varout2(k,idx) = NaN;
        nCorrectAnswers = NaN;
        
    end
    
    if bDebug

        % figure; plot(yi1,PDF1), hold on; 
        % title(sprintf('Target level = %.1f dB, sigma=%.4f',sum_dB_arit([lvl1 lvl2]), sigma))
        % plot(yi2,PDF2,'r')
        % grid on
        % xlim([lvl1-15 lvl1+5]);
        % ylim([0 1])
        % 
        % fprintf('dprime = %.4f, noise normal = %s, noise+test normal = %s\n',dprime(count_row,idx),stats1.bNormal,stats2.bNormal);

        fprintf('noise1 = %.3f (bNormal %.0f), noise2 = %.3f (bNormal %.0f)\n', ...
                                                                    mean( yn1 ), ...
                                                                    Is_normal_distributed(yn1), ...
                                                                    mean( yn2 ), ...
                                                                    Is_normal_distributed(yn2));
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

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inline functions:
function outsig = Il_adjust_tone(insig,fs,lvl)

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(200e-3,fs); ...
               outsigtmp;
               Gen_silence(200e-3,fs)];
           