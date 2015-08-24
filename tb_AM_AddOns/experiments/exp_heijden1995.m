function exp_heijden1995(fs)
% function exp_heijden1995(fs)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/08/2015
% Last update on: 20/08/2015 
% Last use on   : 20/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    fs = 44100; % Hz
end

% Noise generation:
tbuffer = 10;
tn = 500e-3;
fc = 1300;
SPL0 = 60; 
BW1 = 20;
BW2 = 100;

S = Create_sin(fc,tbuffer,fs);
S = setdbspl(S,SPL0);

G1 = AM_random_noise_BW(fc,BW1,SPL0,tbuffer,fs);
G2 = AM_random_noise_BW(fc,BW2,SPL0,tbuffer,fs);

M1 = Multiplied_noise(fc,BW1,SPL0,tbuffer,fs);
M2 = Multiplied_noise(fc,BW2,SPL0,tbuffer,fs);

BP = AM_random_noise(500,800,SPL0-25,tbuffer,fs); % background noise to avoid undesired cues

Sn = S + BP;
G1n = G1 + BP;
G2n = G2 + BP;
M1n = M1 + BP;
M2n = M2 + BP;

% Tone generation:
ft = 2000;
tt = 400e-3;
lenramp = 20;
lensil = (tn - tt)/2;

yt = Create_sin(ft,tt,fs);
yt = setdbspl(yt,SPL0);
yt = Do_cos_ramp(yt,fs,lenramp); 
yt = [Gen_silence(lensil,fs); yt; Gen_silence(lensil,fs)];

model = 'dau1996';
% model = 'jepsen2008';
rampupdown = 20; % ms

opts.model = model;
opts.testlevels      = [];
opts.testlevelsnoise = [60 84];
AboveThres           = [30 75]-SPL0; % taken well above thresholds in Heijden1995, Fig 1. Relative level

testgain = opts.testlevelsnoise - SPL0;

opts.ramp = rampupdown;
opts.siltime = 0;

opts.fmin = fc*0.5; % one octave below masker frequency
opts.fmax = ft*2; % one octave above masker frequency

opts.tmin = 0;
opts.tmax = tn;

Thres = nan(5,length(testgain));

for i = 1:length(testgain)
    opts.gainN = testgain(i);
    opts.Criterion = 1.26;
    opts.AboveThres = AboveThres(i);
    opts.sigmaTimes = 10;
    
    outsS1 = r20150821_update_decision_proc_exp(Sn,yt,opts);
    Thres(1,i) = outsS1.JNDcurrent;
        
    outsG1 = r20150821_update_decision_proc_exp(G1n,yt,opts);
    Thres(2,i) = outsG1.JNDcurrent;

    outsG2 = r20150821_update_decision_proc_exp(G2n,yt,opts);
    Thres(3,i) = outsG2.JNDcurrent;
    
    outsM2 = r20150821_update_decision_proc_exp(M1n,yt,opts);
    Thres(4,i) = outsM2.JNDcurrent;
    
    outsM2 = r20150821_update_decision_proc_exp(M2n,yt,opts);
    Thres(5,i) = outsM2.JNDcurrent;
    
    disp('')
end

if nargout == 0
    
    filename{1} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-S'];
    Wavwrite(Sn,fs,filename{1});
    
    filename{2} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-G1'];
    Wavwrite(G1n,fs,filename{2});
    
    filename{3} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-G2'];
    Wavwrite(G2n,fs,filename{3});
    
    filename{4} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-M1'];
    Wavwrite(M1n,fs,filename{4});
    
    filename{5} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-M2'];
    Wavwrite(M2n,fs,filename{5});
    
    filename{6} = [Get_TUe_paths('outputs') 'exp-heijden-1995-masker-S'];
    Wavwrite(yt,fs,filename{5});
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
