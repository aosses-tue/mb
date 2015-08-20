function y = r20150821_update(f,testtype,model)
% function y = r20150821_update(f,testtype,model)
%
% 1. Description:
%
% 2. Stand-alone example:
%       r20150821_update(1000);
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
    case 3
        dur     = 0.8; % in seconds
        y       = .5*Create_sin(f,dur,fs,0);
        rampup  = 125; % ms
        rampdn  = 125; % ms
        insig    = Do_cos_ramp(y,fs,rampup,rampdn);
        siltime = 0; %200e-3;
        bInsig_is_buffer = 0;
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
sigmaValues = 0.9; % 0.42;%[.3  .5:.1:1]; % [0 0.1 0.5]; % MU

Nsigma      = length(sigmaValues);
JNDcalc     = nan(Nlevels,Nsigma);
CondN       = Nsigma * Nlevels;
bDeterministicTemplate = 0; % if 0, then template is stochastic

for i = 1:Nsigma

    crit            = 70.7;
    
    for j = 1:Nlevels % each testLevel
           
        opts = [];
        opts.f          = f;
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

% ttt = [];
% if bCalculateSigma
%     for j = 1 % for the reference tone
% 
%         lvl = test_level(j);
%         lvl2 = subtract_dB(lvl,0.1);
%         lvl3 = subtract_dB(lvl,0.55);
% 
%         insig1 = setdbspl(y,lvl);
%         insig2 = setdbspl(y,lvl);
%         insig3 = setdbspl(y,lvl3);
%         sigma = 85;
%         
%         outsigsupra         = jepsen2008preproc(gaindb(insig1,5), fs,'resample_intrep'); % 5 dB above the standard level
%         [outsig1, fc, mfc]  = jepsen2008preproc(insig1          , fs,'resample_intrep');
%         outsig2             = jepsen2008preproc(insig2          , fs,'resample_intrep');
%         outsig3             = jepsen2008preproc(insig1+insig3   , fs,'resample_intrep');
% 
%         % idx = 14; 
%         idx = find(fc > 0.5*ft & 2*ft);
%         
%         outsig1_int = []; 
%         outsig4_int = []; 
%         template    = [];
%         for k = 1:50
%             tmp1 = [];
%             tmp4 = [];
%             for i = 1:length(idx)
%                 tmp1 = [tmp1; Add_gaussian_noise(outsig1{idx(i)}(:,1),mu,sigma)];
%                 tmp4 = [tmp4; Add_gaussian_noise(outsigsupra{idx(i)}(:,1), mu, sigma)];
%             end
%             
%             outsig1_int = [outsig1_int tmp1]; 
%             outsig4_int = [outsig4_int tmp4]; 
%         end
%         
%         outsig1_int = mean(outsig1_int,2);
%         outsig4_int = mean(outsig4_int,2);
%         template = [template Get_template_append(outsig1_int,outsig4_int,fs/4)];
%         
%         for k = 1:100
%             
%             outsig1_int = []; 
%             outsig2_int = [];
%             outsig3_int = [];
%             outsign_int = [];
% 
%             for i = 1:length(idx)
%                 tmp1 = Add_gaussian_noise(outsig1{idx(i)}(:,1),mu,sigma);
%                 tmp2 = Add_gaussian_noise(outsig2{idx(i)}(:,1),mu,sigma);
%                 tmp3 = Add_gaussian_noise(outsig3{idx(i)}(:,1),mu,sigma);
%                 tmpn = Add_gaussian_noise(outsig1{idx(i)}(:,1),mu,sigma);
%                 outsig1_int = [outsig1_int; tmp1(:)];
%                 outsig2_int = [outsig2_int; tmp2(:)];
%                 outsig3_int = [outsig3_int; tmp3(:)];
%                 outsign_int = [outsign_int; tmpn(:)];
%             end
%             mue(1) = optimaldetector(outsig1_int,template);
%             mue(2) = optimaldetector(outsig2_int,template);
%             mue(3) = optimaldetector(outsig3_int,template);
% % 
%             [xx maxidx] = max(abs(mue));
%             ttt = [ttt; maxidx mue];
%         end
%         
%         [sum(ttt(:,1)==3) sum( ttt(:,4)>ttt(:,2) ) sum( ttt(:,4)>ttt(:,3) )]
%         
%         crit = abs( max( ttt(:,2),ttt(:,3) )- ttt(:,4) );
%         sum(  crit > sigma )
%         
%     end
% 
% end
% 
% JND = JND(2:end);
% test_level = test_level(2:end);
% 
% figure;
% plot(test_level,JND,'ro','LineWidth',2); grid on; hold on

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
