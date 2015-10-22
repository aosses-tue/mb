function r20151016_update_dau1997_jepsen2008(bParts)
% function r20151016_update_dau1997_jepsen2008(bParts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 13/10/2015
% Last update on: 13/10/2015 
% Last use on   : 13/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;

dir_where  = [Get_TUe_paths('outputs') 'audio-20151006' delim];
dir_where7 = [Get_TUe_paths('outputs') 'audio-20151006-forward' delim];
dir_out    = [Get_TUe_paths('lx_Text') 'lx2015-10-16-decision-CASP' delim 'Figures-new' delim];

%%% Creating stimuli for Fig 2
fname1 = [dir_where 'NBN-fc-5000-Hz-BW-3-Hz-65-dB'];
fname2 = [dir_where 'NBN-fc-5000-Hz-BW-31-Hz-65-dB'];
bCreate = 0;
if bCreate
    masker = exp_dau1997a(fs,3);
    Mkdir(dir_where);
    Wavwrite(masker,fs,fname1);
    
    masker = exp_dau1997a(fs,4);
    Wavwrite(masker,fs,fname2);
end

%%% Creating stimuli for Fig 3
if bCreate
    [masker BWs] = exp_dau1997b(fs,2);
    for i = 1:size(masker,2)
        fname = sprintf('NBN-fc-1000-Hz-BW-%.0f-Hz-65-dB',BWs(i));
        Wavwrite(masker(:,i),fs,fname);
    end
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-10-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-100-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-800-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-1000-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-5000-Hz-65-dB created
end

if nargin == 0
    bParts = [1 0 0 0 0 1];
end

fs = 44100;

bPart1 = bParts(1); % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart2 = bParts(2); % Simultaneous masking, deterministic
bPart3 = bParts(3);
bPart4 = bParts(4); % plotting results of bPart3
bPart5 = bParts(5); % calibration using 800-Hz tone
bPart6 = bParts(6);

opts.nAnalyser      = 101; % 101 = modfilterbank, 103 - jepsen2008

opts.bDecisionMethod = 2; % 2 - cc; 4 - dprime

switch opts.nAnalyser
    case 100
        opts.sigma   = 0.4; % 0.4 for multi-channel, 1.5 for single channel
                             
    case 101
        opts.sigma   = 0.68;    
                                
        opts.modfiltertype = 'dau1997wLP';

    case {103, 104}
        opts.sigma   = 0.615; % 0.615;
    
end
    
opts.var            = opts.sigma.*opts.sigma;
opts.audio.fs       = fs;
if mod(opts.nAnalyser,1)==0
    opts.MethodIntRep   = 1; % 1 - my method; 2 - using casptemplate.m
else
    opts.MethodIntRep   = 2;
end

opts.Reversals4avg  =  6; % 10

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     =  8;

bDebug = 0;
opts.bDebug = bDebug;

close all
count_saved_figures = 1;
if bPart1
    % To do: save template
    bTones = 1;
    bBBN = ~bTones;
    if bTones
        opts.DurRamps   = 125; % additional cosine ramps
        opts.bUseRamp   = 1; % additional cosine ramps
        opts.bUseRampS  = 1; % additional cosine ramps
        erbc2analyse    = freqtoaud([900 1100],'erb'); % 14 for 1000 Hz (approx.)
        opts.filename1 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-60-dB-dur-800-ms.wav'];
        opts.filename2 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-42-dB-dur-800-ms.wav'];
    end
    if bBBN
        opts.DurRamps   = 0; % additional cosine ramps
        opts.bUseRamp   = 0; % additional cosine ramps
        opts.bUseRampS  = 0;
        erbc2analyse    = freqtoaud([100 8000],'erb'); % 14 for 1000 Hz (approx.)
        % opts.filename1 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-60-dB-dur-500-ms.wav'];
        % opts.filename2 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-42-dB-dur-500-ms.wav']; 
        opts.filename1 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-60-dB-dur-500-ms-2.wav'];
        opts.filename2 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-42-dB-dur-500-ms-2.wav']; 
    end
    
    opts = il_get_freqs(erbc2analyse,opts);
    
    opts.Gain4supra = 16; % dB
    opts.audio.fs   = fs;
    
    opts.StepdB     = 10; 
    opts.StepdBmin  = 0.5;
    
    refSPL  = 60; % [20 30 40 50 60 70];
    
    tic
    for i=1:length(refSPL)
        testSPL = refSPL(i)-18; % 42 dB for 60 dB    
        opts.Gain2file1 = refSPL(i) - 60;
        opts.Gain2file2 = refSPL(i) - 60;
        tr_tmp = AMTControl_cl(opts);

        th_JND(i) = sum_dB_arit([refSPL(i) testSPL+tr_tmp.Threshold]) - refSPL(i);
        
        disp('')
    end
    toc
    
    if nargout > 0
        th = th_JND; 
        return;
    end
end

close all
opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart2
    
    opts.DurRamps   = 200; % additional cosine ramps
    opts.Gain4supra =  -3; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 1; % additional cosine ramps
    opts.increment_method = 'modulation-depth';
    
    if strcmp(opts.increment_method,'modulation-depth')
        opts.fmod = 120; % Hz
        opts.dur_test = 1; % s
    end

    opts.Ntimes = 5;
    opts.Nsim = 5;
    opts.filename1 = [dir_where 'NBN-fc-5000-Hz-BW-31-Hz-65-dB.wav'];
    tmp = AMTControl_cl(opts);
    
    if nargout > 0
        th = th_det;
        return;
    end
end

%%%
if bPart3
    
    erbc2analyse    = freqtoaud([4000],'erb'); % 14 for 1000 Hz (approx.) 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-masker-60-dB-dur-10-s.wav'], ... % on-frequency masker
               [dir_where7 'jepsen2008-fig7-2400-Hz-tone-masker-60-dB-dur-10-s.wav'] };   % off-frequency masker
    fnames  = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav'], ...
               [dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-230-ms.wav']}; 
   
	testlevels      = [40 60 70 80; ... % levels for on-freq
                       60 70 80 85];    % levels for off-freq
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.audio.fs   =  fs;
    
    opts.StepdB = 10; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    opts.Silence2noise = 500e-3; % 200 + silence = 700
    opts.increment_method = 'level';
    
    opts.Ntimes = 8;
    opts.Nsim = 1;
    opts.bDebug = 1;
    %%% k = 1: offset-onset of  0 ms
    %%% k = 2: offset-onset of 30 ms
    k = 1;
    opts.filename2 = fnames{k};
    
    for i = 1:length(testlevels) 
        for j = 1
            opts.filename1 = fnamesM{j}; % 1 is on-freq, 2 is off-freq
            opts.Gain2file1 = testlevels(j,i)-60;
            opts.Gain2file2 = opts.Gain2file1;
            tmp = AMTControl_cl(opts);
     
            TTh(j,i) = median(tmp.Threshold)+testlevels(j,i);
            
            disp('')
        end
    end
    
    disp('')
end


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
                  10.00 16.25 25  28.75];    % off-freq, 30-ms, longer interval
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

if bPart5
    
    erbc2analyse    = freqtoaud([800],'erb'); % 14 for 1000 Hz (approx.) 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-800-Hz-tone-masker-60-dB-dur-10-s.wav']}; % on-frequency masker
    fnames  = {[dir_where7 'jepsen2008-fig7-800-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav']}; 
   
	testlevels      = [40 60 70 80]; % levels for on-freq
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.audio.fs   =  fs;
    
    opts.StepdB = 10; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    opts.Silence2noise = 500e-3; % 200 + silence = 700
    opts.increment_method = 'level';
    
    opts.Ntimes = 8;
    opts.Nsim = 1;
    opts.bDebug = 1;
    %%% k = 1: offset-onset of  0 ms
    k = 1;
    opts.filename2 = fnames{k};
    
    for i = 1:length(testlevels) 
        for j = 1
            opts.filename1 = fnamesM{j}; % 1 is on-freq
            opts.Gain2file1 = testlevels(j,i)-60;
            opts.Gain2file2 = opts.Gain2file1;
            tmp = AMTControl_cl(opts);
     
            TTh(j,i) = median(tmp.Threshold)+testlevels(j,i);
            
            disp('')
        end
    end
    
    disp('')
end

if bPart6
    
    testdurs = [10 20 40 70 150]*1e-3*fs;
    fmaskerbuf  = [dir_where 'masker-buf-j2008.wav'];
    fsignals{1} = [dir_where 'sig1-2kHz.wav'];
    fsignals{2} = [dir_where 'sig2-2kHz.wav'];
    fsignals{3} = [dir_where 'sig3-2kHz.wav'];
    fsignals{4} = [dir_where 'sig4-2kHz.wav'];
    fsignals{5} = [dir_where 'sig5-2kHz.wav'];
    
    try
        Wavread(fmaskerbuf{1});
        bCreate = 0;
    catch
        bCreate = 1;
    end
    
    if bCreate
        
        nFig = 4;
        [masker,insig] = exp_jepsen2008(fs,nFig); 
        Wavwrite(masker,fs,fmaskerbuf);
        
        for i = 1:length(testdurs)
            Wavwrite(insig(:,i),fs,fsignals{i});
        end
        
    end
    
    if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 65 dB + 10 = 75 dB
    opts.audio.fs   =  fs;
    
    opts.StepdB    = 10; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 5; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 0;
    opts.increment_method = 'level';
    
    erbc2analyse    = freqtoaud([1000 4000],'erb');
    opts            = il_get_freqs(erbc2analyse,opts);
     
    opts.Ntimes = 8;
    opts.Nsim   = 3;
    opts.bDebug = 1;
    
    opts.filename1 = fmaskerbuf;
        
    for i = 1:length(testdurs) 
        
        opts.filename2 = fsignals{i};
        % opts.Gain2file1 = testlevels(j,i)-60;
        % opts.Gain2file2 = opts.Gain2file1;
        tmp = AMTControl_cl(opts);

        TTh(i) = median(tmp.Threshold)+65;

        disp('')
        
    end
    
    disp('')
    
    
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
