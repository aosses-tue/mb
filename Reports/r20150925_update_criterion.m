function th = r20150925_update_criterion(bParts)
% function th = r20150925_update_criterion(bParts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/09/2015
% Last update on: 02/10/2015 
% Last use on   : 02/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all
outputdir = Get_TUe_paths('outputs');

if nargin == 0
    bParts = [0 0 0 1 1 0 0];
end

fs = 44100;

bPart1 = bParts(1); % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart2 = bParts(2); % Simultaneous masking, deterministic
bPart3 = bParts(3); % Simultaneous masking, stochastic
bPart4 = bParts(4); % Signal integration, deterministic
bPart5 = bParts(5); % Signal integration, stochastic
bPart6 = bParts(6); % Backward masking, deterministic
bPart7 = bParts(7); % Intensity discrimination with BBN

opts.nAnalyser      = 101; % 99 - dau1996a, 99.1, 100 - dau1996, my template estimation

opts.bDecisionMethod = 2; % 2 - cc; 4 - dprime
switch opts.bDecisionMethod
    case 2
        switch opts.nAnalyser
            case 99
                opts.sigma   = 3.25; % tone: 3.25 = band 13-15; 2.7 = band 13-14; 1.9 = band 14; 
            case 100
                opts.sigma   = 0.4; % 1.35; % tone: 3 = band 13-15; % 2.45 = band 13-14; % 1.72 = band 14;
                                     %   BW:                                        1.35 = band 14; (target = 0.83 = -2 dB);
            case 101
                opts.sigma   = 0.68; % 5.8;  % tone: 3.45 = band 14 (all); 1.88 = band 14 (1,2); 1.38 = band 14 (1)
                opts.modfiltertype = 'dau1997wLP';
                
            case 103
                opts.sigma   = 0.615;  %   BW:  21 = band 2-33 (all); 1.95 = band 14; (target = 0.83 = -2 dB);
        end
    case 3
        opts.sigma   = 0.85; % dprime NOT GIVING RELIABLE RESULTS
    case 4
        opts.sigma   = 1.12; % dprime
end
    
opts.var            = opts.sigma.*opts.sigma;
opts.audio.fs       = fs;
opts.MethodIntRep   = 1; % 1 - my method; 2 - using casptemplate.m

opts.Reversals4avg  =  6; % 10

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     =  8;

erbc2analyse        = freqtoaud([1000 1100],'erb'); %freqtoaud([500 2000],'erb'); % 14 for 1000 Hz (approx.)  
opts.fc2plot_idx    = ceil(erbc2analyse(1))-2;
opts.fc2plot_idx2   = floor(erbc2analyse(end))-2;

% dir_where = Get_TUe_paths('outputs');
dir_where = [Get_TUe_paths('outputs') 'audio-20150928' delim];
dir_where2 = [Get_TUe_paths('outputs') 'new' delim 'Experiment' delim];

bDebug = 0;
opts.bDebug = bDebug;

count_saved_figures = 1;
if bPart1
    % To do: save template
    opts.DurRamps   = 150; % additional cosine ramps
    opts.bUseRamp   = 1; % additional cosine ramps
    opts.bUseRampS  = 1; % additional cosine ramps
    opts.Gain4supra =   5; % dB
    opts.audio.fs   =  fs;
    
    opts.StepdB = 2; 
    opts.StepdBmin = 0.2;
    
    opts.filename1 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-60-dB-dur-800-ms.wav'];
    opts.filename2 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-42-dB-dur-800-ms.wav'];
    
    refSPL  = 60; % [20 30 40 50 60 70 80];
    
    for i=1:length(refSPL)
        testSPL = refSPL(i)-18; % 42 dB for 60 dB    
        opts.Gain2file1 = refSPL(i) - 60;
        opts.Gain2file2 = refSPL(i) - 60;
        tr_tmp = AMTControl_cl(opts);

        th_JND(i) = sum_dB_arit([refSPL(i) testSPL+tr_tmp.Threshold]) - refSPL(i);
    end
    
    if nargout > 0
        th = th_JND; % 60 -> 0.6639
        return;
    end
end

close all
opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart2 | bPart3
    
    fmaskerdet{1}  = [dir_where 'masker-det-1.wav'];
    fmaskerdet{2}  = [dir_where 'masker-det-2.wav'];
    fmaskerdet{3}  = [dir_where 'masker-det-3.wav'];
    fmaskerdet{4}  = [dir_where 'masker-det-4.wav'];
    fmaskerdet{5}  = [dir_where 'masker-det-5.wav'];
    idxstart = [0 10 20 50 100]*1e-3*fs;
    fmaskerbuf  = [dir_where 'masker-buf.wav'];
    fsignals{1} = [dir_where 'sig1.wav'];
    fsignals{2} = [dir_where 'sig2.wav'];
    fsignals{3} = [dir_where 'sig3.wav'];
    fsignals{4} = [dir_where 'sig4.wav'];
    fsignals{5} = [dir_where 'sig5.wav'];
    
    try
        Wavread(fmaskerdet{1});
        bCreate = 0;
    catch
        bCreate = 1;
    end
    
    if bCreate
        nFig = 14;
        [masker,insig] = exp_dau1996b(fs,nFig); 

        Wavwrite(masker,fs,fmaskerbuf);
        for i = 1:length(idxstart)
            
            maskerd = [masker(size(insig,1)+1:size(insig,1)+idxstart(i)); masker(1:size(insig,1)-(idxstart(i)))];
            Wavwrite(maskerd,fs,fmaskerdet{i});
            Wavwrite(insig(:,i),fs,fsignals{i});
        end
    end
    
    if bCreate
        il_check_testones(fmaskerdet,fsignals);
    end
    %%%
    
    opts.DurRamps   =  0; % additional cosine ramps
    opts.Gain4supra =  10; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 0; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
end

if bPart2
    opts.Ntimes = 1;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deterministic thresholds:
    for i = 1:length(idxstart)
        opts.filename1 = fmaskerdet{i};
        opts.filename2 = fsignals{i};
        tmp = AMTControl_cl(opts);
        if bDebug
            figure; 
            subplot(2,1,1); plot(tmp.Staircase);  hold on; plot([0 length(tmp.Staircase)],[tmp.Threshold tmp.Threshold],'LineWidth',2)
            subplot(2,1,2); plot(tmp.StaircaseCC')
            close
        end
        th_det(i) = tmp.Threshold;
        disp('')
    end
    
    if nargout > 0
        th = th_det;
        return;
    end
        
    figHandles = findobj('Type','figure'); % 20 generated figures
    figHandles = figHandles(end:-1:1);
    % 1 - noise + tone waveform
    % 2 - noise
    % 3 - templates
    % 4 - staircase
    if bDebug
        h_template = figHandles(3:4:end); % Figure '3 - templates' pooled
    else
        h_template = figHandles;
    end
    plotopts.I_Ylim = [-7 12];
    plotopts.I_Xlim = [0 0.3];
    plotopts.I_Width = 25;
    plotopts.I_Height = 25;
    h2save = Figure2paperfigureT2(h_template,5,1,plotopts);
    xlabel('Time[s]')
    
    Save_all_figures(h2save,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
end

if bPart3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic thresholds:
    close all
    
    opts.Ntimes = 25;
    opts.filename1 = fmaskerbuf;
    Nsims = 4;
    for i = 1:5
        for j = 1:Nsims
            if j == 1
                opts.do_template = 1;
            else
                opts.do_template = 0;
                opts.audio = tmp.audio;
            end
            opts.filename2 = fsignals{i};
            tmp = AMTControl_cl(opts);
            th_sto(j,i) = tmp.Threshold;
            
            if bDebug
                figure; 
                subplot(2,1,1); plot(tmp.Staircase);  hold on; plot([0 length(tmp.Staircase)],[tmp.Threshold tmp.Threshold],'LineWidth',2)
                subplot(2,1,2); plot(tmp.StaircaseCC'); hold on;
                ttmp = tmp.StaircaseCC(3,:)-max(tmp.StaircaseCC(1:2,:));
                plot(ttmp,'r','LineWidth',2)
                disp('')
                close
            end
            
        end
        
    end
    
    if nargout > 0
        th = th_sto;
        return;
    end
    
    figHandles = findobj('Type','figure'); % 20 generated figures
    figHandles = figHandles(end:-1:1);
    
    if bDebug
        h_template = figHandles(3:9:end); % Figure '3 - templates' pooled
    else
        h_template = figHandles;
    end
    plotopts.I_Ylim = [-7 12];
    plotopts.I_Xlim = [0 0.3];
    plotopts.I_Width = 25;
    plotopts.I_Height = 25;
    h2save = Figure2paperfigureT2(h_template,5,1,plotopts);
    xlabel('Time[s]')
    
    Save_all_figures(h2save,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
    
    figure;
    if bPart2
        plot(th_det+75,'ro-'); hold on
    end
    [th_sto_median,errorL,errorU] = Prepare_errorbar_perc(th_sto',25,75);
    % errorbar(median(th_sto,2)+75,'b<--','LineWidth',2); grid on
    errorbar(1:5,th_sto_median+75,errorL,errorU,'b<--','LineWidth',2); grid on
    ha = gca;
    set(gca,'XTick',1:length(th_det))
    set(gca,'XTickLabel',[0 10 20 50 100])
    xlabel('Time onset after masker onset [ms]')
    if bPart2
        legend('Deterministic','Stochastic')
    else
        legend('Stochastic')
    end
    
    Save_all_figures(gcf,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
end

count_saved_figures = 4; warning('temporal')
if bPart4
    
    close all
    nFig = 3;
    [masker,insig] = exp_dau1996b(fs,nFig);
    fmaskerdet{1}  = [dir_where 'masker-det-10.wav'];
    fmaskerdet{2}  = [dir_where 'masker-det-20.wav'];
    fmaskerdet{3}  = [dir_where 'masker-det-30.wav'];
    fmaskerdet{4}  = [dir_where 'masker-det-40.wav'];
    fmaskerdet{5}  = [dir_where 'masker-det-50.wav'];
    dur = [10 20 40 70 150]; % ms
    fmaskerbuf  = [dir_where 'masker-buf0.wav'];
    fsignals{1} = [dir_where 'sig10.wav'];
    fsignals{2} = [dir_where 'sig20.wav'];
    fsignals{3} = [dir_where 'sig30.wav'];
    fsignals{4} = [dir_where 'sig40.wav'];
    fsignals{5} = [dir_where 'sig50.wav'];

    try
        Wavread(fmaskerdet{1});
        bCreate = 0;
    catch
        bCreate = 1;
    end
        
    if bCreate
        Wavwrite(masker,fs,fmaskerbuf);
        for i = 1:length(dur)
            
            maskerd = masker(1:size(insig,1));
            Wavwrite(maskerd,fs,fmaskerdet{i});
            Wavwrite(insig(:,i),fs,fsignals{i});
        end
    end
    
    if bCreate
        il_check_testones(fmaskerdet,fsignals);
    end
    
    opts.DurRamps   = 0; % additional cosine ramps
    opts.Ntimes     = 1;
    opts.bUseRamp   = 1; % additional cosine ramps
    opts.bUseRampS  = 1; % additional cosine ramps
    opts.Gain4supra = 10; % 75 dB + 10 = 85 dB
    opts.audio.fs   = fs;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deterministic thresholds:
    for i = 1:length(dur)
        opts.filename1 = fmaskerdet{i};
        opts.filename2 = fsignals{i};
        
        opts.do_template = 1;
        opts.do_simulation = 1;
        tmp = AMTControl_cl(opts);
        th_det_int(i) = tmp.Threshold;
        
    end
    
    if nargout > 0
        th = th_det_int;
        return;
    end
    
    figHandles = findobj('Type','figure'); % 20 generated figures
    figHandles = figHandles(end:-1:1);
    % 1 - noise + tone waveform
    % 2 - noise
    % 3 - templates
    % 4 - staircase
    if bDebug
        h_template = figHandles(3:4:end); % Figure '3 - templates' pooled
    else
        h_template = figHandles;
    end
    plotopts.I_Ylim = [-3 3.9];
    plotopts.I_Xlim = [0 0.3*12];
    plotopts.I_Width = 25;
    plotopts.I_Height = 25;
    h2save = Figure2paperfigureT(h_template,5,1,plotopts);
    xlabel('Time[s]')
    
    Save_all_figures(h2save,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
    
end

if bPart5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic thresholds:
    close all
    
    opts.Ntimes = 25;
    opts.filename1 = fmaskerbuf;
    Nsims = 4;
    for i = 1:5
        for j = 1:Nsims
            if j == 1
                opts.do_template = 1;
            else
                opts.do_template = 0;
                opts.audio = tmp.audio;
            end
            opts.filename2 = fsignals{i};
            tmp = AMTControl_cl(opts);
            th_sto_int(i,j) = tmp.Threshold;
            
            disp('')
        end
        disp('')
    end
    
    if bDebug
        h_template = figHandles(3:9:end); % Figure '3 - templates' pooled
    else
        h_template = figHandles;
    end
    plotopts.I_Ylim = [-7 12];
    plotopts.I_Xlim = [0 0.3];
    plotopts.I_Width = 25;
    plotopts.I_Height = 25;
    h2save = Figure2paperfigureT2(h_template,5,1,plotopts);
    xlabel('Time[s]')
    
    Save_all_figures(h2save,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
    
    figure;
    if bPart4
        plot(th_det_int+75,'ro-'); hold on
    end
    
    [th_sto_int_median,errorL,errorU] = Prepare_errorbar_perc(th_sto_int',25,75);
    errorbar(1:5,th_sto_int_median+75,errorL,errorU,'b<--','LineWidth',2); grid on
    ha = gca;
    set(ha,'XTick',1:length(th_det_int))
    set(ha,'XTickLabel',dur)
    xlabel('Signal duration [ms]')
    if bPart4
        legend('Deterministic','Stochastic')
    else
        legend('Stochastic')
    end
    
    Save_all_figures(gcf,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
    
end

if bPart6
    
    fmaskerdet{1}  = [dir_where2 'dau1996b_expI_noisemasker-50-ms.wav'];
    
    try
        [x fs] = Wavread(fmaskerdet{1});
    catch
        fmaskertmp  = [dir_where2 'dau1996b_expI_noisemasker'];
        [x fs] = Wavread([fmaskertmp '.wav']);
        x = [Gen_silence(50e-3,fs); x];
        Wavwrite(x,fs,fmaskerdet{1});
        
    end
    % idxstart = [0 10 20 50 100]*1e-3*fs;
    fsignals{1} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-30-ms.wav'];
    fsignals{2} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-35-ms.wav'];
    fsignals{3} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-40-ms.wav'];
    fsignals{4} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-45-ms.wav'];
    fsignals{5} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav'];
    fsignals{6} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-55-ms.wav'];
    fsignals{7} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-60-ms.wav'];
    fsignals{8} = [dir_where2 'dau1996b_expIB0_stim-10ms-76-onset-65-ms.wav'];
    
    tonsets = [-20 -15 -10 -5 0 5 10 15];
    
    t = ( 1:length(x) )/fs;
    figure; plot(t,x); grid on
       
    opts.DurRamps = 0; % additional cosine ramps
    opts.Ntimes = 1;
    opts.Gain4supra =  10; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 0; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    
    opts.filename1 = fmaskerdet{1};
    % Deterministic thresholds:
    for i = 1:length(tonsets)
        opts.filename2 = fsignals{i};
        
        opts.do_template = 1;
        opts.do_simulation = 1;
        tmp = AMTControl_cl(opts);
        th_det_pre(i) = tmp.Threshold;
                
    end
    
    if nargout > 0
        th = th_det_pre;
        return;
    end
    
    figure;
    plot(th_det_pre + 75,'ro-','LineWidth',2); grid on; hold on
    xlabel('Signal onset relative to masker onset [ms]')
    ylabel('Masked threshold [dB]')
    
    set(gca,'XTick',1:length(tonsets));
    set(gca,'XTickLabel',tonsets);
    
end

if bPart7
    % To do: save template
    opts.DurRamps   =   0; % additional cosine ramps
    opts.bUseRamp   =   0; % additional cosine ramps
    opts.bUseRampS  =   0; % additional cosine ramps
    opts.Gain4supra =   5; % dB
    opts.audio.fs   =  fs;
    opts.Ntimes = 1;
    
    opts.StepdB = 2; 
    opts.StepdBmin = 0.2;
    
    opts.filename1 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'jepsen2008-BW-at-60-dB-dur-500-ms.wav'];
    opts.filename2 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'jepsen2008-BW-at-42-dB-dur-500-ms.wav'];
    
    refSPL  = 60; % [20 30 40 50 60 70 80];
    
    for i=1:length(refSPL)
        testSPL = refSPL(i)-18; % 42 dB for 60 dB    
        opts.Gain2file1 = refSPL(i) - 60;
        opts.Gain2file2 = refSPL(i) - 60;
        tr_tmp = AMTControl_cl(opts);

        th_JND(i) = sum_dB_arit([refSPL(i) testSPL+tr_tmp.Threshold]) - refSPL(i);
    end
    
    if nargout > 0
        th = th_JND; % 60 -> 0.6639
        return;
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function il_check_testones(fmaskerdet,fsignals)

N = length(fmaskerdet);
figure;
for i = 1:N
    [m fs] = Wavread(fmaskerdet{i});
    s = Wavread(fsignals{i});
    t = (1:length(m))/fs;
    subplot(N,1,i)
    plot(t,m,t,s,'r'); grid on
end
disp('')
