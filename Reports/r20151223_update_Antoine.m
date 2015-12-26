function y = r20151223_update_Antoine
% function y = r20151223_update_Antoine
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/12/2015
% Last update on: 22/12/2015 
% Last use on   : 22/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bSave   = 1;
bCreate = 0;
hFig    = [];

bComparingICRA_piano = 1;

dir = [Get_TUe_paths('ICRA_Tobias')];
addpath([dir delim 'Tools']);

if bComparingICRA_piano
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sounds to be processed:
    note_test = 'A4';
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];
    
    fname1suffix = 'GRAF28-A4_3';
    fname2suffix = 'JBS51-4544-A4_4';
    fname1   = [dir_where fname1suffix '.wav'];
    fname2   = [dir_where fname2suffix '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dir_out   = [Get_TUe_paths('outputs') 'piano-sounds-for-Antoine' delim];
    Mkdir(dir_out);
    
    [signal1 fs] = Wavread(fname1); 
    [signal2 fs] = Wavread(fname2); 
    
    t1 = ( 1:length(signal1) )/fs;
    t2 = ( 1:length(signal2) )/fs;
    
    [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,10e-3);
    [RMS_se2 t_se2] = rmsdb_sec(signal2,fs,10e-3);
    
    [max_1 idx_1] = max(RMS_se1);
    [max_2 idx_2] = max(RMS_se2);
    delta_dB = max_2 - max_1;
    
    if delta_dB >= 0 % then signal2 is louder
        signal2 = From_dB(-abs(delta_dB+10)) * signal2; % attenuation of the louder signal
    else
        signal1 = From_dB(-abs(delta_dB)) * signal1;
    end
    signal1 = From_dB(-10) * signal1;
    signal2 = From_dB(-10) * signal2;
    
    idx_1 = find(t1 <= t_se1(idx_1),1,'last');
    idx_2 = find(t2 <= t_se2(idx_2),1,'last');
    
    samples_diff = abs(idx_2 - idx_1);
    if idx_1 < idx_2 
    % then idx_2 is after
       [signal2 t2] = Do_alignment(t2,signal2,t2(samples_diff));
    else
    % then idx_2 is after
       [signal1 t1] = Do_alignment(t1,signal1,t1(samples_diff));
    end
    
    [L idxL] = min([length(signal1) length(signal2)]);
    signal1 = Do_truncate(signal1,L);
    signal2 = Do_truncate(signal2,L);
    
    dur_ramp = 150; % ms
    if idxL == 1
    % then signal2 was truncated
        signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp);
    else
    % then signal1 was truncated
        signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp);
    end
    
    clear t1 t2;
    t = (1:L) / fs;
    
    [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,10e-3);
    [RMS_se2 t_se2] = rmsdb_sec(signal2,fs,10e-3);
    
    fc = 20;
    yenv1 = il_get_envelope(signal1,fs,fc);
    yenv2 = il_get_envelope(signal2,fs,fc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    Text_ylabel = 'Amplitude';
    Text_xlabel = 'Time [s]';
    subplot(211)
    plot(t,signal1); hold on, grid on
    plot(t,  yenv1,'k','LineWidth',2);
    plot(t, -yenv1,'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('1','env (approx.)');
    ha = gca;
    
    subplot(212)
    plot(t, signal2,'r'); hold on, grid on
    plot(t,  yenv2,'k','LineWidth',2);
    plot(t, -yenv2,'k','LineWidth',2);
    ylabel(Text_ylabel)
    xlabel(Text_xlabel)
    legend('2','env (approx.)');
    ha(end+1) = gca;
    
    hFig(end+1) = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Copying the input files to the output folder:
    if bSave
        fname_out = sprintf('%s%s-att.wav',dir_out,fname1suffix);
        Wavwrite(signal1,fs,fname_out);
        fname_out = sprintf('%s%s-att.wav',dir_out,fname2suffix);
        Wavwrite(signal2,fs,fname_out);
    end
    
    %%%
    method    = 3; % i = 3 = ERB;
    fname_out1 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname1suffix,method);
    fname_out2 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname2suffix,method);
    
    if bCreate
        noise1   = icra5_noise4piano(signal1,fs,method);
        Wavwrite(noise1,fs,fname_out1);
        noise2   = icra5_noise4piano(signal2,fs,method);
        Wavwrite(noise2,fs,fname_out2);
    else
        noise1  = Wavread(fname_out1);
        noise2  = Wavread(fname_out2);
    end

    % Calculate LTAS
    [ltass1P, fHzP] = calcLTAS(signal1,fs);
    [ltass1N      ] = calcLTAS(noise1,fs);

    % Calculate MTFs
    [mtf1P,cfModHzP] = calcMTF(signal1,fs);
    [mtf1N         ] = calcMTF(noise1,fs);
    
    % Calculate LTAS
    [ltass2P, fHzP] = calcLTAS(signal2,fs);
    [ltass2N      ] = calcLTAS(noise2,fs);

    % Calculate MTFs
    [mtf2P,cfModHzP] = calcMTF(signal2,fs);
    [mtf2N         ] = calcMTF(noise2,fs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    subplot(121);
    semilogx(   fHzP,ltass1P,'-b', ... 
                fHzP,ltass1N,'-k'); hold on, grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({'1','noise'},'location','southwest')
    axis tight;
    XTick = [50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    ha0 = gca;    
    title(name2figname(fname1suffix));
    
    %%%
    subplot(122);
    semilogx(   fHzP,ltass2P,'-r', ... 
                fHzP,ltass2N,'-k'); hold on, grid on;
    xlabel('Frequency [Hz]')
    ylabel('LTASS [dB]')
    legend({'2','noise'},'location','southwest')
    axis tight;
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    ha0(end+1) = gca;
    title(name2figname(fname2suffix));
    
    linkaxes(ha0,'xy');
    ylim([-80 -25])
    xlim([80 8500]);

    hFig(end+1) = gcf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nenv1 = il_get_envelope(noise1,fs,fc);
    nenv2 = il_get_envelope(noise2,fs,fc);

    noise3 = 0.5*noise1+0.5*noise2;
    fname_out3 = sprintf('%snoise-mixed-ICRA-meth-%.0f.wav',dir_out,method);
    
    if bSave
        Wavwrite(noise3,fs,fname_out3);
    end
    nenv3 = il_get_envelope(noise3,fs,fc); % Mix of the two previous noises
    
    figure;
    Text_ylabel = 'Amplitude';
    Text_xlabel = 'Time [s]';
    subplot(311)
    plot(t,noise1); hold on, grid on
    plot(t,  nenv1,'k','LineWidth',2);
    plot(t, -nenv1,'k','LineWidth',2);
    % plot(t_se1, From_dB(RMS_se1,10),'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('1');
    ha(end+1) = gca;
    
    subplot(312)
    plot(t, noise2,'r'); hold on, grid on
    plot(t,  nenv2,'k','LineWidth',2);
    plot(t, -nenv2,'k','LineWidth',2);
    % plot(t_se2, From_dB(RMS_se2,10),'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('2');
    ha(end+1) = gca;
    
    subplot(313)
    plot(t,noise3,'m'); hold on, grid on
    plot(t,  nenv3,'k','LineWidth',2);
    plot(t, -nenv3,'k','LineWidth',2);
    % plot(t_se1, From_dB(RMS_se1,10),'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('mixed');
    ha(end+1) = gca;
    
    linkaxes(ha,'xy');
    xlabel(Text_xlabel)
    ylim([-0.049 0.049])
    xlim([0 max(t)])
    hFig(end+1) = gcf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if bSave
        sig1_SNRm5 = signal1 + From_dB( 5) * noise3;
        sig1_SNRm0 = signal1 + From_dB( 0) * noise3;
        sig1_SNRp5 = signal1 + From_dB(-5) * noise3;
        
        fname_out1 = sprintf('%s%s+mixed-noise-SNR-m05-dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRm5,fs,fname_out1);
        fname_out1 = sprintf('%s%s+mixed-noise-SNR-m00-dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRm0,fs,fname_out1);
        fname_out1 = sprintf('%s%s+mixed-noise-SNR-p05-dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRp5,fs,fname_out1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    subplot(311)
    plot(t, To_dB(yenv2),'r','LineWidth',2); hold on
    plot(t, To_dB(From_dB(5)*nenv3),'k'); grid on
    ha1 = gca;
    title('SNR = -5 dB');
    
    subplot(312)
    plot(t, To_dB(yenv2),'r','LineWidth',2); hold on
    plot(t, To_dB(nenv3),'k'); grid on
    ha1(end+1) = gca;
    title('SNR = 0 dB')
    ylabel('Relative amplitude [dB]')
    
    subplot(313)
    plot(t, To_dB(yenv2),'r','LineWidth',2); hold on
    plot(t, To_dB(From_dB(-5)*nenv3),'k'); grid on
    ha1(end+1) = gca;
    xlabel(Text_xlabel);
    title('SNR = 5 dB')
    
    hFig(end+1) = gcf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MaxdB = max( To_dB(nenv3) ); % Max at 0 dB of SNR
    
    lvlfloor = 20;
    noisefloor = lvlfloor*ones(size(t));
    env2plot = max( To_dB(nenv3), MaxdB-lvlfloor); 
    
    figure;
    plot(t, To_dB(yenv2),'r','LineWidth',2); hold on; grid on;
    plot(t, env2plot,'k','LineWidth',2);
    plot(t, MaxdB-noisefloor,'k--')
    ha1(end+1) = gca;
    
    linkaxes(ha1,'xy');
    xlim([0 max(t)])
    ylim([-65 -30])
    xlabel(Text_xlabel);
    ylabel('Relative amplitude [dB]')
    hFig(end+1)=gcf;
    disp('')
end

if bSave
    for i = 1:length(hFig)
        fname = sprintf('%sICRA-fig-%.0f',dir_out,i);
        Save_figure_as(hFig(i),fname,'epsc');
        Save_figure_as(hFig(i),fname,'fig');
    end
end

rmpath([dir delim 'Tools'])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);
