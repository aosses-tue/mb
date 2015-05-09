function r20150501_update
% function r20150501_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 29/04/2015
% Last update on: 08/05/2015 % Update this date manually
% Last use on   : 08/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bPart1 = 0;
bPart2 = 0; % Determining the noise deviation: Weber's approach
bPart3 = 1; % Determining the noise deviation: Weber's approach
bPart4 = 1;

% Common parameters:
close all
model = 'dau1996a';
fs  = 44100;

Delta = 1/fs;
t = ( 1:44100 )/fs;
S = ones(length(t),1);
sum( Delta*(S.*S) );


if bPart1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f       = 3150;
dur     = 10e-3;
lvl1    = 87; % Suprathreshold
lvl2    = 80;
lvl3    = lvl2+1;

lvlnoise= 77;

ynoiseb = AM_random_noise(20,5000,lvlnoise,200e-3,fs);
ynoiseb = setdbspl(ynoiseb,lvlnoise);

ynoise = [  Gen_silence(100e-3,fs); ...
            ynoiseb;
            Gen_silence(300e-3,fs)];

y       = Create_sin(f,dur,fs,0);

rampup = 5; % ms
rampdn = 5; % ms
y1  = Do_cos_ramp(y,fs,rampup,rampdn);
y1b = setdbspl(y1,lvl1);
y2b = setdbspl(y1,lvl2);
y3b = setdbspl(y1,lvl3);

y1 = [  Gen_silence(200e-3,fs); ...
        y1b;
        Gen_silence(390e-3,fs)];
y2 = [  Gen_silence(200e-3,fs); ...
        y2b;
        Gen_silence(390e-3,fs)];
y3 = [  Gen_silence(200e-3,fs); ...
        y3b;
        Gen_silence(390e-3,fs)];


SMT = ynoise + y1;  % Supra-threshold
SM  = ynoise;       % Masker alone
SMTc= ynoise + y2;  % Current signal
SMTc2= ynoise + y3;  % Current signal

% sound(y1 ,fs)
% sound(SMT,fs)

setup.bAddNoise = 0;
setup.fc = f;
% setup.LvNoise   = 30;
[out  fc t] = Get_internal_representations_deterministic([SMT SM SMTc SMTc2],fs,model,setup);

RMT = out(:,1);
RM  = out(:,2);
RMTc= out(:,3);
RMTc2 = out(:,4);

setup = [];
setup.fs = fs;
Template = Get_template(RM,RMT,setup);

% cc1 = Get_decision_criterion(RMTc -RM,Template,'cross-correlation');
% cc2 = Get_decision_criterion(RMTc2-RM,Template,'cross-correlation');
cc1 = Correlation(RMTc -RM,Template,'coefficient');
cc2 = Correlation(RMTc2-RM,Template,'coefficient');

% cc1 = Get_decision_criterion(template,RMTc ,'cross-correlation',opts);
% cc2 = Get_decision_criterion(template,RMTc2,'cross-correlation',opts);

% label1 = mean(outsig1(end-N+1:end));
% label2 = mean(outsig2(end-N+1:end));
% Delta = label2 - label1;

figure;
subplot(3,1,1)
plot(t,RM); grid on
xlabel('Time [s]')

subplot(3,1,2)
plot(t,RMT), grid on %, t,outsig2); grid on
xlabel('Time [s]')
legend

subplot(3,1,3)
plot(t,RMT-RM), grid on %, t,outsig2); grid on
xlabel('Time [s]')
legend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% teoJND      = [1.44 1.13 .95 .82  .80]; % dB
% teoLevels   = [20   40   60   80  100]; % Level test tone
teoJND      = [1.6   .95 1.]; % dB
teoLevels   = [20   60  100]; % Level test tone


f       = 1000; %3150;
dur     = 1; % in seconds

y       = .5*Create_sin(f,dur,fs,0);
rampup  = 5; % ms
rampdn  = 5; % ms
ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);

for idx = 3; % 1:length(teoLevels) % idx     = 3; % related to 60 dB
    lvl1    = teoLevels(idx); 
    lvl2    = subtract_dB( lvl1, teoJND(idx) );
    lvl3    = teoLevels(idx)+3;
    
    y1b     = setdbspl(ytmp,lvl1);
    y2b     = setdbspl(ytmp,lvl2);
    y3b     = setdbspl(ytmp,lvl3);

    y1 = [  Gen_silence(200e-3,fs); ...
            y1b;
            Gen_silence(200e-3,fs)];
    y2 = [  Gen_silence(200e-3,fs); ...
            y2b;
            Gen_silence(200e-3,fs)];
    y3 = [  Gen_silence(200e-3,fs); ...
            y3b;
            Gen_silence(200e-3,fs)];
        
    SMT     = y3;
    SMTc    = y1;  % Current signal
    SMTc2   = y2;  % Current signal

    % sound(y1 ,fs)
    % sound(y2 ,fs)

    switch f
        case 1000
            idxband = 14;
        case 3150
            idxband = 23;
        otherwise
            error('put manually the fc idx manually');
    end
    
    setup.bAddNoise = 0;
    setup.sigma     = 10; % set here the deviation of the noise
    setup.fs        = fs;
    % [out  fc t] = Get_internal_representations_deterministic([SMT SMTc SMTc2],fs,model,setup);
    % RMT = out(:,1);
    % RMTc= out(:,2);
    % RMTc2 = out(:,3);
    
    [RMT    fc t] = Get_internal_representations(SMT  ,fs,model,setup);
    [RMTc   fc t] = Get_internal_representations(SMTc ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations(SMTc2,fs,model,setup);
    RMT     = RMT(:,idxband);
    RMTc    = RMTc(:,idxband);
    RMTc2   = RMTc2(:,idxband);

    N = length(RMT);
    mu = 0;
    yn1 = normrnd(mu,setup.sigma,N,1);
    yn2 = normrnd(mu,setup.sigma,N,1);
    yn3 = normrnd(mu,setup.sigma,N,1);
     
    setup = [];
    setup.fs = fs;
    template = Get_template(RMTc,RMT,setup); % RMTc is here noise alone
%     cc(idx) = Get_decision_criterion(RMTc, RMTc2 ,'cross-correlation',setup);
    cc(idx) = Get_decision_criterion(template, RMTc2 ,'cross-correlation',setup);
    
    cct(1)  = Get_decision_criterion(template    , template     ,'cross-correlation',setup);
    cct2(1) = Get_decision_criterion(template    , template     ,'cross-correlation-non-normalised',setup);
    cct(2)  = Get_decision_criterion(template+yn1, template+yn2 ,'cross-correlation',setup);
    cct2(2) = Get_decision_criterion(template+yn1, template+yn2 ,'cross-correlation-non-normalised',setup);
    
    difference(idx) = 1/fs*sum(RMTc2 - RMTc); % current difference

    disp('')
end

fprintf('Correlation of: %.2f\n',cc);
fprintf('Difference of: %.2f\n',difference);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes respect bPart2:
% 
%   - Estimate threshold for N times
%   - No suprathreshold signal
%   - No template used

stepJND     = 0.2;
testJND     = [1.6 1.1 .7]; %[2:-stepJND:.8 0]; % 0.8:stepJND:1.0; % dB
% testJND     = [1.5 1.1 1 .9 .75 0]; 
testLevels  = 60*ones(size(testJND)); % Level test tone

sigmaValues = [0  20  30  35];
sigmaTimes  = [1 100 100 100];

f       = 1000; 
dur     = 4; % in seconds

y       = .5*Create_sin(f,dur,fs,0);
rampup  = 5; % ms
rampdn  = 5; % ms
ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);

durtot  = dur+400e-3;
% idx_compare = 1:round(durtot*fs);
idx_compare = round(1*fs):round(2*fs);

for idx = 1:length(testJND) 

    lvl1    = testLevels(idx);
    SMTc    = Il_create_tone(ytmp,fs,lvl1);
    
    lvl2    = subtract_dB( lvl1, testJND(idx) );
    SMTc2   = Il_create_tone(ytmp,fs,lvl2);
        
    % sound(SMTc ,fs)
    % sound(SMTc2,fs)

    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    
    [RMTc   fc t] = Get_internal_representations_deterministic(SMTc      ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);
    N = length(RMTc);
    
    count_row = 1;

    for idx_s = 1:length(sigmaValues)
        mu = 0;
        sigma = sigmaValues(idx_s);
        
        for idx_times = 1:sigmaTimes(idx_s)
            yn1 = normrnd(mu,sigma,N,1);
            yn2 = normrnd(mu,sigma,N,1);

            RMTc_n  = RMTc  + yn1;
            RMTc2_n = RMTc2 + yn2;
                
            if count_row == 1
                figure;
                plot( 1/fs*(RMTc2_n(idx_compare) - RMTc_n(idx_compare)) ); grid on
            end
            difference(count_row,idx) = 1/fs*sum(RMTc2_n(idx_compare) - RMTc_n(idx_compare)); % current difference
            count_row = count_row + 1;
            
        end
    end
    
end

M0dB = mean( difference(:,end) );
% difference = difference/M0dB;

idCum = cumsum(sigmaTimes);

diff1 = difference(1,:);
diff2 = difference(idCum(2-1)+1:idCum(2),:);
diff3 = difference(idCum(3-1)+1:idCum(3),:);
diff4 = difference(idCum(4-1)+1:idCum(4),:);

crit = 70.7; % percentage correct

P1  = percentile([diff1; diff1],crit);
P2  = percentile(diff2,crit);
P3  = percentile(diff3,crit);
P4  = percentile(diff4,crit);

bigP = [P1;P2;P3;P4]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same as 3, but now it looks for a threshold in levels other than 60 dB SPL:

if ~bPart3
    error('Set bPart3 first');
end

idx_n = 3;
Criterion   = []; %P3(idx_n);
txt = sprintf('Criterion = interp1(testJND,P%.0f,0.95);',idx_n);
eval(txt);

sigma       = sigmaValues(idx_n);
times       = sigmaTimes(idx_n);
% testJND     = [1 .82      .6  0]; %[0.95 .82]; % 0.8:stepJND:1.0; % dB
testLevels  = 100*ones(size(testJND));

for idx = 1:length(testJND) 

    lvl1    = testLevels(idx);
    SMTc    = Il_create_tone(ytmp,fs,lvl1);
    
    lvl2    = subtract_dB( lvl1, testJND(idx) );
    SMTc2   = Il_create_tone(ytmp,fs,lvl2);
        
    % sound(SMTc ,fs)
    % sound(SMTc2,fs)

    setup.bAddNoise = 0;
    setup.fs        = fs;
    setup.fc        = f;
    
    [RMTc   fc t] = Get_internal_representations_deterministic(SMTc ,fs,model,setup);
    [RMTc2  fc t] = Get_internal_representations_deterministic(SMTc+SMTc2,fs,model,setup);
    N = length(RMTc);
    
    count_row = 1;

    for idx_times = 1:times
        
        yn1 = normrnd(mu,sigma,N,1);
        yn2 = normrnd(mu,sigma,N,1);

        RMTc_n  = RMTc  + yn1;
        RMTc2_n = RMTc2 + yn2;

        diff80(count_row,idx) = 1/fs*sum(RMTc2_n(idx_compare) - RMTc_n(idx_compare)); % current difference
        count_row = count_row + 1;

    end
    
end

M0dB = mean( diff80(:,end) );
%diff80 = diff80/M0dB;
P3new  = percentile(diff80,crit);

amount = sum(diff80<=Criterion);
tmp = interp1(P3new,testJND,Criterion);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


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
           