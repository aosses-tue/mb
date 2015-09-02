function [hFig data] = exp_heijden1995(fs,bParts)
% function [hFig data] = exp_heijden1995(fs,bParts)
%
% 1. Description:
%
% 2. Stand-alone example:
%       % To generate Figures 1 (no error bars) and 2 of Heijden1995. The 
%       % figure handles are returned in hFig.
%       fs      = 44100; 
%       bParts  = [0 1 0]; 
%       hFig = exp_heijden1995(fs,bParts);
% 
%       fs = 44100;
%       bParts  = [0 0 1]; 
%       hFig = exp_heijden1995(fs,bParts);
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/08/2015
% Last update on: 01/09/2015 
% Last use on   : 01/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    fs = 44100; % Hz
end

if nargin < 2
    bParts = [0 1 0];
end

if nargout == 2
    bParts(2) = 1; 
end

bGenSignals         = bParts(1);
bGenOriginalPlots   = bParts(2);
bPlotSimulation     = bParts(3);

hFig = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bGenSignals
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bGenOriginalPlots
    
    FontSize = 14;
    
    testLevels = 60:6:84;
    Tr = [  14 23 36 54 65.5; ...     % Sine tone masker
            12 18 27 40 52.5; ...     % 100-Hz GN
            12 16 23.5 36.5 49; ...   % 20-Hz GN
            12 17 23 35 45; ...       % 100-Hz MN
            11.5  14 21 29 40];       % 20-Hz MN
	
	labels = {'T','100-Hz GN','20-Hz GN','100-Hz MN','20-Hz MN'};
    
	figure;
    plot(testLevels,Tr(1,:),'bo-','LineWidth',1), hold on
    plot(testLevels,Tr(2,:),'rs--','LineWidth',2)
    plot(testLevels,Tr(3,:),'ks-.','LineWidth',1)
    plot(testLevels,Tr(4,:),'r>-.','LineWidth',2)
    plot(testLevels,Tr(5,:),'k>-','LineWidth',1)
    legend(labels,'Location','NorthWest')
    grid on
    xlabel('Masker level [dB SPL]','FontSize',FontSize)
    ylabel('Threshold [dB SPL]','FontSize',FontSize)
    xlim([59 86])
    hFig(end+1) = gcf;
    set(gca,'FontSize',FontSize);
    
    figure;
    plot(testLevels,Tr(1,:)-Tr(2,:),'rs--','LineWidth',2), hold on
    plot(testLevels,Tr(1,:)-Tr(3,:),'ks-.','LineWidth',1)
    plot(testLevels,Tr(1,:)-Tr(4,:),'r>-.','LineWidth',2)
    plot(testLevels,Tr(1,:)-Tr(5,:),'k>-','LineWidth',1)
    legend(labels{2:end},'Location','NorthWest')
    grid on
    xlabel('Masker level [dB SPL]','FontSize',FontSize)
    ylabel('Masking release [dB]','FontSize',FontSize)
    xlim([59 86])
    hFig(end+1) = gcf;
    set(gca,'FontSize',FontSize);
    
    data.Thresholds = Tr;
    data.signaltypes = labels;
    data.testLevels = testLevels;
end

if bPlotSimulation
    
    FontSize = 14;
    
    testLevels = 60:12:84;
    Tr = [  -27.5 -24.5 -16; ...   % Sine tone masker
            -34.8 -26.1 -22.8; ...   % 100-Hz GN
            -41.5  -32  -26; ...   % 20-Hz GN
            -41.9  -28  -25; ...     % 100-Hz MN
            -43.5  -31.5 -33.5];      % 20-Hz MN
	Tr(:,1) = Tr(:,1) + 60;
    Tr(:,2) = Tr(:,2) + 72;
    Tr(:,3) = Tr(:,3) + 84;
    
	labels = {'T','100-Hz GN','20-Hz GN','100-Hz MN','20-Hz MN'};
    
	figure;
    plot(testLevels,Tr(1,:),'bo-','LineWidth',1), hold on
    plot(testLevels,Tr(2,:),'rs--','LineWidth',2)
    plot(testLevels,Tr(3,:),'ks-.','LineWidth',1)
    plot(testLevels,Tr(4,:),'r>-.','LineWidth',2)
    plot(testLevels,Tr(5,:),'k>-','LineWidth',1)
    legend(labels,'Location','NorthWest')
    grid on
    xlabel('Masker level [dB SPL]','FontSize',FontSize)
    ylabel('Threshold [dB SPL]','FontSize',FontSize)
    xlim([59 86])
    hFig(end+1) = gcf;
    set(gca,'FontSize',FontSize);
    
    figure;
    plot(testLevels,Tr(1,:)-Tr(2,:),'rs--','LineWidth',2), hold on
    plot(testLevels,Tr(1,:)-Tr(3,:),'ks-.','LineWidth',1)
    plot(testLevels,Tr(1,:)-Tr(4,:),'r>-.','LineWidth',2)
    plot(testLevels,Tr(1,:)-Tr(5,:),'k>-','LineWidth',1)
    legend(labels{2:end},'Location','NorthWest')
    grid on
    xlabel('Masker level [dB SPL]','FontSize',FontSize)
    ylabel('Masking release [dB]','FontSize',FontSize)
    xlim([59 86])
    hFig(end+1) = gcf;
    set(gca,'FontSize',FontSize);
    
    data.Thresholds = Tr;
    data.signaltypes = labels;
    data.testLevels = testLevels;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
