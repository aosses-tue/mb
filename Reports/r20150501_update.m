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
% Last update on: 29/04/2015 % Update this date manually
% Last use on   : 29/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bPart1 = 0;
bPart2 = 1; % Determining the noise deviation: Weber's approach

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

teoJND      = [1.44 1.13 .95 .82  .80]; % dB
teoLevels   = [20   40   60   80  100]; % Level test tone

f       = 1000; %3150;
dur     = 1; % in seconds
idx     = 1; % related to 60 dB
lvl1    = teoLevels(idx); 
lvl2    = lvl1 + teoJND(idx);
lvl3    = lvl1 + 25;

y       = .5*Create_sin(f,dur,fs,0);

rampup  = 5; % ms
rampdn  = 5; % ms
ytmp    = Do_cos_ramp(y,fs,rampup,rampdn);
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

SMT     = y3;  % Supra-threshold
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
setup.bAddNoise = 1;
setup.sigma     = 100; % set here the deviation of the noise

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


setup = [];
setup.fs = fs;
template = Get_template(RMT*0,RMT,setup);
cc1 = Get_decision_criterion(template,RMTc ,'cross-correlation',setup);
cc2 = Get_decision_criterion(template,RMTc2,'cross-correlation',setup);

difference = 1/fs*sum(RMTc2 - RMTc);

figure;
plot(RMTc), hold on
plot(RMTc2, 'r')

figure;
plot(SMTc-SMTc2)

% % outsig2 = outsig2(:,idx);
% N = length(RS1)/10;
% % N = length(outsig1);
% 
% % label1 = mean(outsig1(end-N+1:end));
% % label2 = mean(outsig2(end-N+1:end));
% % Delta = label2 - label1;

figure;
subplot(2,1,1)
plot(t,RMTc,t,RMTc2); grid on
legend('lvl1','lvl1+JND')
xlabel('Time [s]')
title(sprintf('Channel tuned to fc=%.0f [Hz], lvl = %.2f [dB]',fc(idxband),lvl1))
subplot(2,1,2)
plot(t,RMTc-RMTc2), grid on %, t,outsig2); grid on
xlabel('Time [s]')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
