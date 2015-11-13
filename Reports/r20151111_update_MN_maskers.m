function r20151111_update_MN_maskers
% function r20151111_update_MN_maskers
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/11/2015
% Last update on: 12/11/2015 
% Last use on   : 12/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

fs = 44100;

dir_where7 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim];
dir_out    = [Get_TUe_paths('lx_Text') 'lx2015-10-16-decision-CASP' delim 'Figures-new' delim];

if nargin == 0
    %         1 2 3 4 5 6 7 8
    bParts = [1 0 1 0 0 0 1 0];
    % CC:     0 1 0 - - 1 0 -
    %                   0       % if 2000-Hz  
end

fs = 44100;

bPart3 = bParts(3); % Forward masking: on and off-frequency
bPart4 = bParts(4); % plotting results of bPart3

opts.nAnalyser      = 103; % 101 = modfilterbank, 103 - jepsen2008

opts.bDecisionMethod = 5; % 2 - cc; 4 - dprime; 5 - cc updated

switch opts.nAnalyser
    case 101
        opts.sigma   = 0.8; % 0.68;    
                                
        opts.modfiltertype = 'dau1997wLP';

    case {103, 104}
        opts.sigma   = 0.44; %0.88; % 0.59; % 0.615 = 17.86 dB; % 0.45 if gain_after_drnl = 13 dB
    
end
    
opts.var            = opts.sigma.*opts.sigma;
opts.audio.fs       = fs;
opts.MethodIntRep   = 1; % 1 - my method; 2 - using casptemplate.m

opts.Reversals4avg  = 6; % 10

opts.do_template    = 1;
opts.do_simulation  = 1;
opts.Nreversals     = 10;
opts.nDown          = 2;
opts.experiment_type = 'AFC';
% opts.experiment_type = 'constant';

bDebug = 0;
opts.bDebug = bDebug;

close all
count_saved_figures = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart3
    
    erbc2analyse    = freqtoaud([3500 5000],'erb'); % 14 for 1000 Hz (approx.)  % 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'sine-1300-Hz-60-dB.wav'], ...
               [dir_where7 'mult-noise-1300-Hz-BW-20-Hz-60-dB.wav'], ... 
               [dir_where7 'mult-noise-1300-Hz-BW-100-Hz-60-dB.wav'], ...
               [dir_where7 'randomnoise-fc-1300_BW-20_fmod-0_Mdept-0_SPL-60.wav'], ...
               [dir_where7 'randomnoise-fc-1300_BW-100_fmod-0_Mdept-0_SPL-60.wav']};
    fnames  = {[dir_where7 'sine-2000-Hz-ramps-of-20-ms-60-dB.wav']}; 
   
	testlevels = 80;
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.Level_start = 10;
    opts.audio.fs   =  fs;
    
    opts.StepdB     = 8; 
    opts.StepdBmin  = 2;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 0;
    % opts.Silence2noise = 500e-3; % 200 + silence = 700. If sil = 400e-3, then signal dur = 300e-3 -> simultaneous masking
    opts.increment_method = 'level';
    
    opts.Ntimes = 8;
    opts.Nsim   = 4;
    opts.bDebug = 0;
    
    for k = 2:5
        opts.filename2 = fnames{1};

        for i = 1:length(testlevels) 
            for j = 1% :2
                
                opts.filename1 = fnamesM{k}; 
                opts.Gain2file1 = testlevels(i)-60;
                opts.Gain2file2 = opts.Gain2file1;
                tmp = AMTControl_cl(opts);

                TTh(k) = median(tmp.Threshold)+testlevels(i);
                
                disp('')
                
            end
        end
    end
    
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart4
    
    hFig4 = []; % handles for the figures generated in this part
    
    %% 4.a Results of intensity discrimination task
    % raw1 - intensity discrimination with 1-kHz tones
    % raw2 - intensity discrimination with BBN (100 Hz - 8 kHz)
    testlevels = [20 30 40 50 60 70];
    minL = 20; maxL = 70;
    minY = 0.3; maxY = 1.5;
    
    raw1_d1997 = [1.3480    1.0299    0.8509    0.7214    0.6281    0.5465]; % sigma = 0.68
    raw2_d1997 = [0.5941    0.5778    0.7214    0.8233    0.8279    0.8279]; % sigma = 0.68
    raw1_j2008 = [0.6109    0.5168    0.4752    0.4752    0.5168    0.5314]; % sigma = 0.615
    raw2_j2008 = [0.6109    0.6639    0.6281    0.7017    0.8279    1.0873]; % sigma = 0.615
    % raw1_j2008 = [0.7623 0.5204]; % sigma = 0.59 on 30/10/2015, SPL = 20, 60
    % raw2_j2008 = [0.7623 1.0023]; % sigma = 0.59; SPL = 60, 70
    
    figure;
    subplot(1,2,1)
    plot(testlevels, raw1_d1997,'rs-','LineWidth',2); hold on
	plot(testlevels, raw1_j2008,'ko-','LineWidth',2); grid on
	xlabel('Standard level [dB]'); xlim([minL maxL])
    ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
    title('A. test 1-kHz tones');
    set(gca,'XTick',testlevels);
    set(gca,'YTick',[minY:0.1:maxY])
    
    subplot(1,2,2)
    plot( testlevels, raw2_d1997,'rs-','LineWidth',2); hold on
    plot( testlevels, raw2_j2008,'ko-','LineWidth',2); grid on
    xlabel('Standard level [dB]'); xlim([minL maxL])
    ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
    title('B. test BBNs');
    set(gca,'XTick',testlevels);
    set(gca,'YTick',[minY:0.1:maxY])
    
    legend('PEMO', 'CASP','Location','NorthWest');
    hFig4(end+1) = gcf;
        
    %%%
    sigmas      = [0.5 0.6 0.615 0.63 0.7 0.85];
    LineColor  = {'bs-';'m>-.';'ko-';'r--';'g>-';'r<-'};
    LineWidth  = [1.5 1 2 1 1 2];
    
    raw3_j2008 = [0.4752    0.4015    0.3796    0.3796    0.4015    0.4248; ... % sigma = 0.5
                  0.5778    0.4887    0.4369    0.4369    0.4621    0.4752; ... % sigma = 0.6
                  0.6109    0.5168    0.4752    0.4752    0.5168    0.5314; ... % sigma = 0.615
                  0.6109    0.5314    0.4752    0.4887    0.5285    0.5465; ... % sigma = 0.63
                  0.6639    0.5619    0.5314    0.5314    0.5619    0.5941; ... % sigma = 0.7
                  0.8279    0.7017    0.6281    0.6639    0.7017    0.7416];    % sigma = 0.85
 
    raw4_j2008 = [0.5026    0.5314    0.5026    0.5619    0.6639    0.8745; ...
                  0.5941    0.6458    0.6281    0.6826    0.8279    1.0583; ...
                  0.6109    0.6639    0.6281    0.7017    0.8279    1.0873; ...
                  0.6281    0.6658    0.6458    0.7214    0.8509    1.1171; ...
                  0.7017    0.7416    0.7214    0.7836    0.9492    1.2441; ...
                  0.8509    0.8988    0.8745    0.9492    1.1476    1.4990];
    
	figure;
    for i = 1:length(LineColor)
        subplot(1,2,1)
        plot(testlevels, raw3_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
        if i==1;
            title('C. test 1-kHz tones');
            xlabel('Standard level [dB]')
            ylabel('\Delta at threshold [dB]')
            xlim([minL maxL])
            ylim([minY maxY])
            set(gca,'XTick',testlevels);
            set(gca,'YTick',[minY:0.1:maxY]);
        end
        subplot(1,2,2)
        plot(testlevels, raw4_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
        if i==1;
            title('D. test BBNs');
            xlabel('Standard level [dB]')
            ylabel('\Delta at threshold [dB]')
            xlim([minL maxL])
            ylim([minY maxY])
            set(gca,'XTick',testlevels);
            set(gca,'YTick',[minY:0.1:maxY]);
        end
        Leg{i} = ['\sigma = ' num2str(sigmas(i))];
    end
    legend(Leg,'Location','NorthWest')
    hFig4(end+1) = gcf;
    
    %% 4.b Results of forward-masking experiment with tones
    testlevels = [  40 60 70 80; ...
                    60 70 80 85];
	minL = min(min(testlevels));
    maxL = max(max(testlevels));
    
	% Obtained on 13-oct-2015:
    raw1_d1997 = [17.5000   32.5000   45.0000   59.0000; ...	% on-freq,  0-ms, longer interval
                  17.5000   32.5000   41.2500   51.2500];       % on-freq, 30-ms, longer interval

	raw2_d1997 = [11.2500   22.5000   33.7500   38.7500; ...	% off-freq,   0-ms, longer interval 
                   6.2500    8.7500   15.0000   17.5000];       % off-freq,  30-ms, longer interval     

    raw1_j2008 = [18.75 24.25 28.125 33.75; ... % 	on-freq,  0-ms, longer interval
                  20.00 26.25 28.750 31.25];    %	on-freq, 30-ms, longer interval];
    % 20.0000   26.2500   28.7500   31.2500	on-freq, 30-ms, longer interval
    % 20.0000   26.2500   28.7500   31.3750	on-freq, 30-ms, longer interval (run 2)

    raw2_j2008 = [21.25 27.5 35.0 45.875; ... %	off-freq,  0-ms, longer interval
                  10.00 16.25 25  28.75];     % off-freq, 30-ms, longer interval
    % 60	70	80	85
    % 21.2500   27.5000   35.0000   46.2500	off-freq,  0-ms, longer interval
    % 21.2500   27.5000   35.0000   45.5000	off-freq,  0-ms, longer interval (run 2)
    
    % ref1_d1997;
    % ref2_d1997;
    % ref1_j2008;
    % ref2_j2008;
    
    figure;
    subplot(1,2,1)
    plot(testlevels(1,:), raw1_d1997(1,:),'r>-','LineWidth',2); hold on
	plot(testlevels(1,:), raw1_d1997(2,:),'ro--');
    plot(testlevels(1,:), raw1_j2008(1,:),'k>-');
	plot(testlevels(1,:), raw1_j2008(2,:),'ko--'); grid on
    xlabel('Masker level [dB SPL]'); xlim([minL maxL])
    ylabel('Threshold [dB SPL]')
    title('A. On-freq.');
    
    subplot(1,2,2)
    plot( testlevels(2,:), raw2_d1997(1,:),'r>-','LineWidth',2); hold on
    plot( testlevels(2,:), raw2_d1997(2,:),'ro--')
    plot( testlevels(2,:), raw2_j2008(1,:),'k>-')
	plot( testlevels(2,:), raw2_j2008(2,:),'ko--'); grid on
    xlabel('Masker level [dB SPL]'); xlim([minL maxL])
    ylabel('Threshold [dB SPL]')
    title('B. Off-freq.');
    
    legend('PEMO,  0-ms', ...
           'PEMO, 30-ms', ...
           'CASP,  0-ms', ...
           'CASP, 30-ms','Location','NorthWest');
    hFig4(end+1) = gcf;
       
	Save_all_figures(hFig4,dir_out,count_saved_figures);
    count_saved_figures = count_saved_figures + length(hFig4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function opts = il_get_freqs(erbc2analyse, opts)

if length(erbc2analyse) == 1
    opts.fc2plot_idx    = round(erbc2analyse(1))-2;
    opts.fc2plot_idx2 = opts.fc2plot_idx;
else
    opts.fc2plot_idx    = ceil(erbc2analyse(1))-2;
    opts.fc2plot_idx2   = floor(erbc2analyse(end))-2;
end