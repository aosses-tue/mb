function r20150925_update_criterion
% function r20150925_update_criterion
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
% Last update on: 22/09/2015 
% Last use on   : 22/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all
outputdir = Get_TUe_paths('outputs');

fs = 44100;

bPart1 = 1; % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart2 = 1; % Deterministic
bPart3 = 1; % Stochastic
bPart4 = 1; % Deterministic
bPart5 = 1; % Stochastic

opts.bDecisionMethod = 2;
opts.sigma        = 1.7; % 0.95; 
opts.audio.fs     = fs;
opts.nAnalyser    = 100; % 99 - dau1996a, 99.1, 100 - dau1996, my template estimation
opts.MethodIntRep = 1; % 1 - my method; 2 - using casptemplate.m

opts.Reversals4avg = 10;

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     = 12;

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
    refSPL  = 60;
    testSPL = 42;
    tr_tmp = AMTControl_cl(opts);
    
    th_JND = sum_dB_arit([refSPL testSPL+tr_tmp.Threshold]) - refSPL;
    
end

close all
opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart2 | bPart3
    
    fmaskerdet{1}  = [Get_TUe_paths('outputs') 'masker-det-1.wav'];
    fmaskerdet{2}  = [Get_TUe_paths('outputs') 'masker-det-2.wav'];
    fmaskerdet{3}  = [Get_TUe_paths('outputs') 'masker-det-3.wav'];
    fmaskerdet{4}  = [Get_TUe_paths('outputs') 'masker-det-4.wav'];
    fmaskerdet{5}  = [Get_TUe_paths('outputs') 'masker-det-5.wav'];
    idxstart = [0 10 20 50 100]*1e-3*fs;
    fmaskerbuf  = [Get_TUe_paths('outputs') 'masker-buf.wav'];
    fsignals{1} = [Get_TUe_paths('outputs') 'sig1.wav'];
    fsignals{2} = [Get_TUe_paths('outputs') 'sig2.wav'];
    fsignals{3} = [Get_TUe_paths('outputs') 'sig3.wav'];
    fsignals{4} = [Get_TUe_paths('outputs') 'sig4.wav'];
    fsignals{5} = [Get_TUe_paths('outputs') 'sig5.wav'];
    
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
        th_det(i) = tmp.Threshold;
    end
    figHandles = findobj('Type','figure'); % 20 generated figures
    figHandles = figHandles(end:-1:1);
    % 1 - noise + tone waveform
    % 2 - noise
    % 3 - templates
    % 4 - staircase
    h_template = figHandles(3:4:end); % Figure '3 - templates' pooled
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
    
    opts.Ntimes = 100;
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
            th_sto(i,j) = tmp.Threshold;
            
            disp('')
        end
        disp('')
    end
    
    h_template = figHandles(3:9:end); % Figure '3 - templates' pooled
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
    plot(median(th_sto,2)+75,'b<--','LineWidth',2); grid on
    ha = gca;
    set(gca,'XTick',1:length(th_det))
    set(gca,'XTickLabel',[0 10 20 50 100])
    xlabel('Time onset after masker onset [ms]')
    if bPart2
        legend('Deterministic','Stochastic')
    else
        legend('Stochastic')
    end
    
end

if bPart4
    
    close all
    nFig = 3;
    [masker,insig] = exp_dau1996b(fs,nFig);
    fmaskerdet{1}  = [Get_TUe_paths('outputs') 'masker-det-10.wav'];
    fmaskerdet{2}  = [Get_TUe_paths('outputs') 'masker-det-20.wav'];
    fmaskerdet{3}  = [Get_TUe_paths('outputs') 'masker-det-30.wav'];
    fmaskerdet{4}  = [Get_TUe_paths('outputs') 'masker-det-40.wav'];
    fmaskerdet{5}  = [Get_TUe_paths('outputs') 'masker-det-50.wav'];
    dur = [10 20 40 70 150]; % ms
    fmaskerbuf  = [Get_TUe_paths('outputs') 'masker-buf0.wav'];
    fsignals{1} = [Get_TUe_paths('outputs') 'sig10.wav'];
    fsignals{2} = [Get_TUe_paths('outputs') 'sig20.wav'];
    fsignals{3} = [Get_TUe_paths('outputs') 'sig30.wav'];
    fsignals{4} = [Get_TUe_paths('outputs') 'sig40.wav'];
    fsignals{5} = [Get_TUe_paths('outputs') 'sig50.wav'];

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
    
    opts.DurRamps = 0; % additional cosine ramps
    opts.Ntimes = 1;
    opts.bUseRamp   = 1; % additional cosine ramps
    opts.bUseRampS  = 1; % additional cosine ramps
    opts.Gain4supra =   5; % dB
    opts.audio.fs   =  fs;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deterministic thresholds:
    for i = 1:length(dur)
        opts.filename1 = fmaskerdet{i};
        opts.filename2 = fsignals{i};
        
        opts.do_template = 1;
        opts.do_simulation = 1;
        tmp = AMTControl_cl(opts);
        th_det_int(i) = tmp.Threshold;
        
        disp('')
    end
    
    figHandles = findobj('Type','figure'); % 20 generated figures
    figHandles = figHandles(end:-1:1);
    % 1 - noise + tone waveform
    % 2 - noise
    % 3 - templates
    % 4 - staircase
    h_template = figHandles(3:4:end); % Figure '3 - templates' pooled
    plotopts.I_Ylim = [-7 12];
    plotopts.I_Xlim = [0 0.3];
    plotopts.I_Width = 25;
    plotopts.I_Height = 25;
    h2save = Figure2paperfigureT2(h_template,5,1,plotopts);
    xlabel('Time[s]')
    
    Save_all_figures(h2save,outputdir,count_saved_figures);
    count_saved_figures = count_saved_figures + 1;
    
end

if bPart5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic thresholds:
    close all
    
    opts.Ntimes = 100;
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
    
    h_template = figHandles(3:9:end); % Figure '3 - templates' pooled
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
        plot(th_det_int+75,'ro-'); hold on
    end
    plot(median(th_sto_int,2)+75,'b<--','LineWidth',2); grid on
    ha = gca;
    set(gca,'XTick',1:length(th_det_int))
    set(gca,'XTickLabel',dur)
    xlabel('Signal duration [ms]')
    if bPart2
        legend('Deterministic','Stochastic')
    else
        legend('Stochastic')
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
