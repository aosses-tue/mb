function r20151002_update_dau
% function r20151002_update_dau
%
% 1. Description:
%       This script is a copy of r20150320_update_dauM.m, but re-run in Sept
%       2015, where some parts of the code were optimised (and updated).
% 
% 2. Stand-alone example:
%       r20151002_update_dau;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: r20150320_update, r20150320_updateN
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/09/2015
% Last update on: 28/09/2015 
% Last use on   : 29/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

bDiary = 0;
Diary(mfilename,bDiary);

bPart1          = 0;
bPart2          = 0; % Simulation
bPart3          = 0; % d-prime theory
bPart4          = 1;

bPart2recalculate = 1;
if bPart1
    bCreateAudio= 0; % run this first
end

opts.calc_method = 2;

bPlot = 1;

h               = [];

tmp = Get_TUe_paths('outputs');
pathaudiosrc  = [tmp 'new' delim '48000' delim];
pathaudiodst1 = [tmp 'new' delim '44100' delim];% 'D:\Output\r20141031_update_dau_et_al20150318\';
pathaudiodst2 = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\'; % APEX results
pathaudiodst3 = [tmp 'new' delim 'Experiment' delim];
output_dir = [Get_TUe_paths('lx_Text') 'lx2015-03-20-update' delim 'Figures' delim];
output_figures = [Get_TUe_paths('lx_Text') 'lx2015-09-28-decision' delim 'Figures-new' delim];

test_onset_ref  = [-20:5:10]*1e-3;

x = [];
offset = 0.2;

if bPart1
    % figure;
    for i = 1:10
        fname   = [pathaudiosrc 'dau1996b_expIB0_stim-10ms-76-' num2str(i) '.wav'];
        fname2  = sprintf('%sdau1996b_expIB0_stim-10ms-76-onset-%.0f-ms',pathaudiodst1,test_onset_ref(i)*1000+50);
        [tmp fs]= Wavread(fname);
        x(:,i) = tmp;
        t = (0:length(tmp)-1)/fs;

        % plot(t*1000,tmp+offset*(i-1)), hold on

        y = resample(tmp,44100,fs);
        if bCreateAudio
            Wavwrite(y,44100,fname2);
        end
    end

    N = size(y,1);
    fname = [pathaudiodst1 'dau1996b_expI_noisemasker.wav']; % existing already
    [tmp fs] = Wavread(fname);
    M = size(tmp,1);
    x(:,i+1) = [tmp; tmp(1:N-M-1)];
else
    % fill in filenames
    % generates 50ms of silence prior to the noise and truncates to noise length
end

dirResults = [pathaudiodst2 'Results' delim];
dirStimuli = [pathaudiodst3 delim]; % same as [pathaudiodst2 'Stimuli' delim];

file1   = [dirResults 'dau1996b-IIC-3AFC-multi-1-of-2-results-3.apr.xml']; 
file2   = [dirResults 'dau1996b-IIC-3AFC-multi-1-of-2-results-4.apr.xml'];
file3   = [dirResults 'dau1996b-IIC-3AFC-multi-2-of-2-results-1.apr.xml']; 
file4   = [dirResults 'dau1996b-IIC-3AFC-multi-2-of-2-results-2.apr.xml'];

stim1   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-30-ms.wav']; % -20
stim2   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-35-ms.wav']; % -15
stim3   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-40-ms.wav']; % -10
stim4   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-45-ms.wav']; % -5
stim5   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav']; % 0
stim6   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-55-ms.wav']; % 5
stim7   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-60-ms.wav']; % 10
stimnoise = [dirStimuli 'dau1996b_expI_noisemasker.wav'];

filenames{1} = stim1;
filenames{2} = stim2;
filenames{3} = stim3;
filenames{4} = stim4;
filenames{5} = stim5;
filenames{6} = stim6;
filenames{7} = stim7;
filenames{8} = stimnoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2
    
    if bPart2recalculate
        opts.bExpIC0 = 1; % Backward masking

        opts.method = 'dau1996';
        opts.pathaudio = pathaudiodst3;
        opts.filenames = filenames;
        r20141031_update_dau_et_alN(opts);

        opts.method = 'dau1996a';
        opts.pathaudio = pathaudiodst3;
        opts.filenames = filenames;
        r20141031_update_dau_et_alN(opts);

    else

        error('Not updated yet')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Backward masking:
        % path = 'D:\Output\r20141031_update_dau_et_al20150318\';
        % f1 = [path 'mue-2015-03-18-at-19h-59m-40s.mat']; % no limitation
        % f2 = [path 'mue-2015-03-18-at-20h-11m-24s.mat']; % limitation
        % 
        % load(f1); % loads mue
        % 
        % dB_SPL          = [10 16:10:76]; % Get this from processing
        % test_onset_ref  = [-100:40:-20,-15:5:15]*1e-3; % Get this from processing
        % 
        % criterion_corr = 8.5;% options.criterion_corr;
        % 
        % Thres   = demo_dau1996b_decision(mue(3:end,:),dB_SPL,criterion_corr);
        % 
        % load(f2); % loads mue
        % ThresO  = demo_dau1996b_decision(mue(3:end,:),dB_SPL,criterion_corr);
        % 
        % figure;
        % plot(test_onset_ref(3:end)*1000,Thres,'rx-'), grid on, hold on
        % plot(test_onset_ref(3:end)*1000,ThresO,'bo-','MarkerFaceColor','b');
        % 
        % errorbar(stim_onset,MeanV,StdV,'k-')
        % legend('no limit','limit = 10','Measured Thresholds','Location','NorthWest')
        % xlabel('Signal onset relative to masker onset [ms]')
        % ylabel('Masked thresholds [dB]')
        % 
        % h(end+1) = gcf;
    end
end

if bPart3 % to prove that error function implementation is the same that 
    % equation A8 in Dau1996a
    muvar = 0:0.2:4; 
    prob2 = 1 - (erfc((((muvar)*0.707) - 0    ) * sqrt(2)/2) / 2);
    prob3 = 1 - (erfc((((muvar)*0.765) - 0.423) * sqrt(2)/2) / 2);
    prob4 = 1 - (erfc((((muvar)*0.810) - 0.668) * sqrt(2)/2) / 2);
    figure; 
    plot(   muvar,100*prob2, ...
            muvar,100*prob3, ...
            muvar,100*prob4); grid on; hold on
    
    hold on;
    mu = [0.00012 0.55653 0.83865];
    sigma = [1.41415 1.30384 1.25471];
    for i = 1:3
        z(i,:) = normcdf(  (muvar - mu(i))/sigma(i)  );
    end
    % figure;
    plot(muvar,100*z,'o','LineWidth',2);
    xlabel('Ratio \mu/\sigma or ''d-prime''')
    ylabel('Percent correct [%]')
    
    legend('2-AFC','3-AFC','4-AFC')
end

if bPart4
    
    choice = 'cc';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results intensity discrimination task (results obtained on 02/10/2015):
    %       th_JND obtained from r20150925_update_decision
    
    switch choice
        case 'reference'
            figtitle = 'dau1996';
            
            lvl1ref = [20 40 60 80 100];
            ref1 = [1.44 1.13 0.95 0.82 0.80]; % dB
            
            lvl2plot = lvl1ref;
            th2plot = ref1;
                        
        case 'cc'
            figtitle = 'Crit 1';
            %    To recalculate: % resnew = r20150925_update_criterion([1 0 0 0 0 0 0]);
            lvltest = [20 30 40 50 60 70 80];
            res = [1.1171  1.0299  1.0023  0.9701  0.9492  0.9492  0.9808]; % Method 2, std = 1.75 (thres = std)
            % res=[1.1476  1.0583  1.0023  0.9754  0.9492  0.9544  0.9969]; % Method 2, std = 2.45 (thres = std) x 2-Ch: 13-14
            lvl2plot = lvltest;
            th2plot = res;
            
        case 'd-prime'
            figtitle = 'Crit 2';
            % Method 4, std = 1.12
            lvltest = [20 30 40 50 60 70 80];
            res = [ 1.5149    1.1790    1.0583    0.9969    0.9544    0.9492    0.9237];
            
            lvl2plot = lvltest;
            th2plot = res;
    end
    
    figure; 
    plot(lvl2plot, th2plot,'ko--', 'LineWidth', 2); hold on; grid on
    xlim([20 80])
    xlabel('Test level [dB]')
    ylabel('\Delta L [dB]')
    title(figtitle)
   
    href = gcf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dau1996b, Fig15
    
    switch choice
        case 'reference'
            tonsetref = [0    5    15 35 55 75 95]; % ms
            refdet    = [66.1 67.2 72 74 76 77 78.5]; % at 5 and 35 ms response were interpolated
            refsto    = [79.2 79   79.2 79.5 79 80.5 80]; 
            
            t2plot = tonsetref;
            th2plot(1,:) = refdet;
            th2plot(2,:) = refsto;
            
        case 'cc'
            
            tonsettest = [0 10 20 50 100]; % ms
            
            t2plot = tonsettest;
            
            th2plot = [-11 -5 -11 -12 -14]+75; % resnew = r20150925_update_criterion([0 1 0 0 0 0 0]);
            
            storaw1 = [ -1.5 -2.0  2.0  0.5  3.0; ...
                        -2.0  1.0  1.0  1.0  2.0; ...
                        -2.5 -0.5  4.0  1.0  1.5; ...
                        -2.0  1.0  3.5  1.0  3.0]; % resnew = r20150925_update_criterion([0 0 1 0 0 0 0]);
            th2plot(2,:) = median(storaw1)+75;
                                 
        case 'd-prime'
            error('Calculate')
            
            teststoraw2 = [ 0.5000   -1.0000    1.0000    4.0000; ...
                            4.5000    7.5000    4.5000    4.5000; ...
                            10.0000   13.0000   10.0000   11.5000; ...
                            14.0000   11.5000   13.5000   14.0000; ...
                            16.0000   11.0000   18.0000   15.0000];
            teststo(2,:) = median(teststoraw2')+75;
    end
    
    %%% Plot references:
    figure;
    plot(t2plot, th2plot(1,:),'ko--','MarkerFaceColor','k'); hold on; grid on
    plot(t2plot, th2plot(2,:),'k>--','MarkerFaceColor','k'); 
    legend('deterministic','stochastic')
    xlabel('Signal onset relative to masker onset [ms]')
    ylabel('Masked threshold [dB SPL]')
    title(figtitle)
    href(end+1)=gcf;
	%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dau1996b, Fig3
    switch choice
        case 'reference'
            tdur = [10 20 40 80 160 320 400];
            ref = [71 69 67 64 63.5 64 65];
            
            t2plot = tdur;
            th2plot = ref;
            
        case 'cc'
            tdur = [10 20 40 70 150];
            th2plot = [-18   -13   -25   -27   -23]+75; % resnew = r20150925_update_criterion([0 0 0 1 0 0 0]);
            t2plot = tdur;
            
    end
    
    figure;
    plot(t2plot, th2plot,'ko--','MarkerFaceColor','k'); hold on; grid on
	xlabel('Signal duration [ms]')
    ylabel('Masked threshold [dB SPL]')
    title(figtitle)
    href(end+1) = gcf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dau1996b
    switch choice
        case 'reference'
            tonsetref = [-100 -60 -30 -20 -15 -12.5 -10  -5 -2.5 0 5 10 15 20];
            ref    = [  18 17.5   17  17  18  19    24  48  70 79 75 66 75 74];
            
            t2plot = tonsetref;
            th2plot = ref;
            
        case 'cc'
            tonsets = [-20 -15 -10 -5 0 5 10 15];
            % resnew = r20150925_update_criterion([0 0 0 0 0 1]);
            
            t2plot = tonsets;
            th2plot = [-64   -64   -62   -46   -19   -18   -13   -16]+75;
            
    end
    
    figure;
    plot(t2plot, th2plot,'ko--','MarkerFaceColor','k'); hold on; grid on
	xlabel('Signal onset relative to masker onset [ms]')
    ylabel('Masked threshold [dB SPL]')
    href(end+1)=gcf;
    title(figtitle)
    
    plotopts.FontSize = 14;
    plotopts.I_Width = 20;
    plotopts.I_Height = 8;
    plotopts.I_KeepYTicks = 1;
    href2 = Figure2paperfigureT(href,1,4,plotopts);
    
    switch choice
        case 'reference'
            Save_figure_as(href2,[output_figures 'figs-dau1996b'],'epsc');
        case 'cc'
            Save_figure_as(href2,[output_figures 'figs-crit1'],'epsc');
    end
    
end
Save_all_figures(h,output_dir);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
