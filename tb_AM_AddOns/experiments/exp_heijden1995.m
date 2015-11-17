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
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Data simulated on: XXXXXX
%     testLevels = 60:12:84;
%     % Tr = [  -27.5 -24.5 -16; ...   % Sine tone masker
%     %         -34.8 -26.1 -22.8; ...   % 100-Hz GN
%     %         -41.5  -32  -26; ...   % 20-Hz GN
%     %         -41.9  -28  -25; ...     % 100-Hz MN
%     %         -43.5  -31.5 -33.5];      % 20-Hz MN
%     Tr = [  -56.0  -53.75 -46.5; ...   % Sine tone masker
%             -48.75 -42.75 -34.0; ...   % 100-Hz GN
%             -50.5  -48.75 -34.5; ...   % 20-Hz GN
%             NaN NaN NaN; ...     % 100-Hz MN
%             NaN NaN NaN];      % 20-Hz MN
%     
% 	  Tr(:,1) = Tr(:,1) + 60;
%     Tr(:,2) = Tr(:,2) + 72;
%     Tr(:,3) = Tr(:,3) + 84;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data simulated on: 11/11/2015
    
    bUse1 = 1; % Jepsen2008
    if bUse1
        % testLevels = 60:5:85;
        % Data simulated on: 11/11/2015 (supra at +10dB)
        % Tr =  [32.0000   37.5000   44.5000   49.5000   60.0000   67.0000; ... % Sine tone: fairly good (higher shift at 60 dB)
        %        28.0000   34.5000   41.0000   45.5000   55.0000   62.0000; ... % 100-Hz GN % More difficult than in reality (at 85 dB should be 52 dB) 
        %        30.0000   37.0000   44.5000   51.5000   66.0000   71.0000; ... %  20-Hz GN % Much more difficult than in reality (at 85 dB should be 48 dB)
        %        23.5000   33.5000   35.5000   42.0000   49.5000   55.5000; ... % 100-Hz MN % More difficult than in reality (at 85 dB should be 45 dB) 
        %        27.0000   35.5000   44.0000   52.5000   61.5000   70.0000]; ... %  20-Hz MN % Much more difficult than in reality (at 85 dB should be 40 dB)
        % 
        % TTh25 = [31.0000   36.5000   43.0000   48.0000   59.5000   63.5000; ...
        %          26.5000   33.0000   38.0000   44.0000   53.0000   60.5000; ...
        %          29.5000   33.5000   43.5000   47.0000   63.0000   71.0000; ...
        %          22.5000   32.0000   32.5000   41.0000   48.5000   53.5000; ...
        %          26.0000   34.0000   42.0000   50.0000   59.0000   69.5000];
        % 
        % TTh75 = [32.0000   38.5000   45.0000   51.5000   61.0000   67.5000; ...
        %          28.0000   36.5000   42.0000   47.5000   58.0000   63.0000; ...
        %          31.0000   38.0000   46.5000   54.0000   68.0000   71.5000; ...
        %          26.0000   34.0000   37.0000   44.5000   50.5000   58.5000; ...
        %          29.0000   37.5000   47.0000   54.0000   63.5000   70.5000];
        
        % Data simulated on: 17/11/2015 (supra at -10dB, between 1500-4000)
        testLevels = 60:6:84;
        Tr    = [  28 44   61; ...
                   26.0000   38.0000   55.5000; ...
                   28.0000   47.0000   64.5000; ...
                   21.0000   32.5000   47.5000; ...
                   27.5000   38.5000   65.0000];
        TTh25 = [  25 44   59; ...
                   26.0000   37.0000   54.5000; ...
                   27.0000   44.0000   63.0000; ...
                   20.5000   31.5000   47.0000; ...
                   26.5000   38.0000   64.0000];
        TTh75 = [   28 44.5 62; ...
                   27.0000   38.5000   56.5000; ...
                   28.0000   48.5000   65.5000; ...
                   21.5000   33.5000   49.0000; ...
                   28.0000   39.5000   66.5000];

    else
        testLevels = 60:5:85;
        Tr = [23.5000   27.0000   33.5000   38.5000   40.5000   47.0000; ...
              25.5000   31.0000   35.5000   42.0000   48.0000   52.5000; ...
              29.0000   31.5000   39.5000   43.5000   52.0000   58.0000; ...
              22.5000   28.5000   34.5000   39.0000   46.5000   50.5000; ...
              28.0000   34.0000   42.0000   46.5000   51.0000   56.5000];

        TTh25 = [ 22.5000   26.5000   31.5000   37.5000   36.5000   47.0000; ...
                  24.0000   30.0000   34.0000   39.0000   47.5000   49.5000; ...
                  25.5000   30.0000   38.5000   42.5000   50.5000   55.5000; ...
                  19.0000   26.0000   32.0000   38.0000   45.0000   49.0000; ...
                  27.0000   31.0000   41.0000   44.5000   48.0000   55.5000];
 
        TTh75 =  [24.0000   27.0000   34.5000   39.0000   43.5000   48.0000; ...
                  27.0000   31.5000   38.0000   43.0000   48.0000   55.0000; ...
                  29.5000   33.5000   40.0000   44.5000   53.5000   58.0000; ...
                  24.5000   29.0000   36.0000   40.0000   47.5000   51.0000; ...
                  28.5000   35.5000   42.5000   47.5000   53.0000   57.0000];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
