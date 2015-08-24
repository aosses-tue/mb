function y = r20150821_update(f,testtype,model)
% function y = r20150821_update(f,testtype,model)
%
% 1. Description:
%
% 2. Stand-alone example:
%       r20150821_update(1000,3,'jepsen2008'); % Weber's law as in Jepsen 2008, 1-kHz tones
%       r20150821_update(1000,2,'jepsen2008'); % Weber's law as in Jepsen 2008, BBN
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 17/08/2015
% Last update on: 17/08/2015 
% Last use on   : 17/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bPart1 = 0; % actual modelling part
bPart2 = 1; % Plot results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params:
if nargin < 3
    % model = 'dau1996a';
    model = 'jepsen2008';
end

if nargin < 2
    % testtype = 1; % sine tones as in Dau1996a
    % testtype = 2; % narrow-band as in Jepsen2008 (used to get variance)
    testtype = 3; % sine tones as in Jepsen2008
end

if nargin < 1
    f       = 5000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart1
% stepJND     = [ .04; ...
%                 .04; ...
%                 .005; ...
%                 .04; ...
%                 .04];  

%                                        %  dB    dau1996     dau1997     dau1997multi
% testJNDi    = [  26     1.00 1.40; ... %  26    1.1304      1.1369
%                  40     0.90 1.30; ... %  40    1.0750      1.0788   
%                  60     1.00 1.05; ... %  60    1.0040      1.0049
%                  80     0.60 1.00; ... %  80    0.9311      0.9342
%                 100     0.60 1.00];    % 100    0.8666      0.8657
%                                        %        var=0.8    var=0.42                

stepJND     = [ .04; ...
                .04; ...
                .04; ...
                .04; ...
                .04; ...
                .04];  
%                                        %  dB    dau1996     dau1997     dau1997multi
% testJNDi    = [  26     0.40 0.80; ... %  26    1.1304      1.1369
%                  40     0.40 0.80; ... %  40    1.0750      1.0788   
%                  60     0.50 0.90; ... %  60    1.0040      1.0049
%                  80     0.70 1.10; ... %  80    0.9311      0.9342
%                 100     0.90 1.30];    % 100    0.8666      0.8657
%                                        %        var=0.8    var=0.42                             

testJNDi    = [  20     0.40 0.90; ...
                 30     0.40 0.90; ...
                 40     0.40 0.90; ... %  40    1.0750      1.0788  
                 50     0.40 0.90; ...
                 60     0.50 1.00; ... %  60    1.0040      1.0049
                 70     0.60 1.10; ... %  80    0.9311      0.9342
                 80     0.70 1.20];    % 100    0.8666      0.8657
                                       %        var=0.8    var=0.42                             


tmpJND = testJNDi(1,2):stepJND(1):testJNDi(1,3);    
% NJND        = length(tmpJND);
% testLevels  = [26 40 60 80 100]; 
testLevels  =[30 40 50 60 70]; % [30 40 60 70 80]; 
Nlevels     = length(testLevels);
Ntimes      = 2; %10; % 1 calculations for each condition
CondCounter = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating test tone:
fs      = 44100;

switch testtype
    case 1
        dur     = 0.5; % in seconds
        y       = .5*Create_sin(f,dur,fs,0);
        rampup  = 5; % ms
        rampdn  = 5; % ms
        insig    = Do_cos_ramp(y,fs,rampup,rampdn);
        siltime = 0; %200e-3;
        bInsig_is_buffer = 0;
        
        fmin = f; % to use single-channel model
        fmax = f; % to use single-channel model
    case 2
        dur     = 10; % buffer in seconds
        % BW      = 100;
        % y       = AM_random_noise_BW(f,BW,60,dur,fs); % in the experiment the level is set
        y       = AM_random_noise(100,10000,60,dur,fs); % in the experiment the level is set
        rampup  = 50; % ms
        rampdn  = 50; % ms
        insig   = y;
        bInsig_is_buffer = 1;
        % insig    = Do_cos_ramp(y,fs,rampup,rampdn);
        siltime = 0;
        
        fmin    = 100; % to use multi-channel model, all fcs above 100 Hz
        fmax    = 8000; % to use multi-channel model, all fcs below 8000 Hz
    case 3
        dur     = 0.8; % in seconds
        y       = .5*Create_sin(f,dur,fs,0);
        rampup  = 125; % ms
        rampdn  = 125; % ms
        insig    = Do_cos_ramp(y,fs,rampup,rampdn);
        siltime = 0; %200e-3;
        bInsig_is_buffer = 0;
        
        fmin = 0.5*f; % to use multi-channel model, down to one octave below fc
        fmax = 2*f; % to use multi-channel model, up to one octave above fc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Criterion = 1.26; % D-prime for 3-AFC 

% tmin = 0; tmax = 4;
% sigmaValues = 2;
% model = 'dau1996';
switch testtype
    case 1
        tmin = 0; tmax = 0.5; % 4; % 2 3
    case 2
        tmin = 0; tmax = 500e-3;
    case 3
        tmin = 0; tmax = 800e-3;
end
sigmaValues = 0.95; % 0.42;%[.3  .5:.1:1]; % [0 0.1 0.5]; % MU

Nsigma      = length(sigmaValues);
JNDcalc     = nan(Nlevels,Nsigma);
CondN       = Nsigma * Nlevels;
bDeterministicTemplate = 0; % if 0, then template is stochastic

for i = 1:Nsigma

    crit            = 70.7;
    
    for j = 1:Nlevels % each testLevel
           
        opts = [];
        opts.f          = f;
        opts.fmin       = fmin;
        opts.fmax       = fmax;
        opts.Criterion  = Criterion;
        opts.testLevel  = testLevels(j);
        fprintf('   - variance = %.2f: lvl = %.0f dB \n',sigmaValues(i),testLevels(j));
        
        idxLvl = find( testJNDi(:,i)==testLevels(j) );
        opts.testJND    = testJNDi(idxLvl,2):stepJND:testJNDi(idxLvl,3);
        
        opts.insig      = insig;
        opts.crit       = crit;
        opts.sigma      = sigmaValues(i);
        opts.sigmaTimes = Ntimes;
        opts.bDebug     = 1;
        opts.model      = model;
        opts.tmin = tmin;
        opts.tmax = tmax;
        opts.fs         = fs;
        opts.bDeterministicTemplate = bDeterministicTemplate;
        opts.siltime = siltime;
        opts.bInsig_is_buffer = bInsig_is_buffer;
        outs2           = r20150821_update_decision_proc(opts);
        
        JNDcalc(j,i)    = outs2.JNDcurrent;
        fprintf('Completed %.2f%% (%.0f out of %.0f conditions)\n',CondCounter/CondN*100,CondCounter,CondN);
        CondCounter     = CondCounter + 1;
    end
end

var2latex([testLevels' JNDcalc]);
end % end bPart1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2
    lvl_j2008 = 20:10:70;
    data_j2008_BBN = [0.49 0.55 0.67 0.62 0.83 0.97];
    data_j2008_1kHz = [0.65 0.65 0.52 0.52 0.68 0.68];
                  % 20  30      40      50      60      70
    data_own_BBN = [NaN 0.46981 0.45913 0.50751 0.61146 0.76693; ... 0.85   % NaN 0.47046 0.46591 0.50224 0.6161 0.72225; ... 0.85 old
                    NaN 0.49807 0.48094 0.54841 0.66558 0.76845; ... 0.90   % NaN 0.4977  0.48393 0.54102 0.66287 0.76268; ... 0.90, old
                    NaN 0.52576 0.50963 0.56151 0.692   0.81209]; % 0.95    % NaN 0.52348 0.51517 0.55842 0.6987  0.79719]; % 0.95, old
	data_own_1kHz= [NaN 0.74964 0.69323 0.71756 0.75882 0.79411; ... % 0.85
                    NaN 0.79373 0.73414 0.75991 0.80429 0.84142; ... 0.90
                    NaN 0.83824 0.77472 0.80215 0.84836 0.8878]; % 0.95
    
	compBBN = repmat(data_j2008_BBN,3,1) - data_own_BBN;
    MBBN = mean(compBBN(:,2:end),2);
    comp1kHz = repmat(data_j2008_1kHz,3,1) - data_own_1kHz;
    M1kHz = mean(comp1kHz(:,2:end),2);
    TotDiff = MBBN + abs(M1kHz);
    
    figure;
    plot(lvl_j2008,data_j2008_BBN,lvl_j2008,data_own_BBN)
end

%         outsigsupra         = jepsen2008preproc(gaindb(insig1,5), fs,'resample_intrep'); % 5 dB above the standard level

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
