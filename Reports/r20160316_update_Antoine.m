function r20160316_update_Antoine
% function r20160316_update_Antoine
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

bSave   = 0;
bCreate = 1;
hFig    = [];

bPart1  = 0;
bPart2  = 0; 
bPart3  = 0; % Analysis of the filter bank
bPart4  = 1; % L338. Generating/explaining the pedestal noise

dir = [Get_TUe_paths('ICRA_Tobias')];
addpath([dir delim 'Tools']);

note_test = 'A4';
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];

dir_out   = [Get_TUe_paths('outputs') 'piano-sounds-for-Antoine' delim];
Mkdir(dir_out);
dir_out_figs = [dir_out 'Figures' delim];
Mkdir(dir_out_figs);
    
if bPart1 || bPart2

    piano_1 = 'GRAF28';
    piano_2 = 'JBS51-4544';
    attenuate_by = 10; % Because these signals were measured very closed to the soundboard (they are very loud)
    
    % 1. Sounds to be processed:
    fname1suffix = [piano_1 '-' note_test '_3'];
    fname2suffix = [piano_2 '-' note_test '_4'];
    fname1   = [dir_where fname1suffix '.wav'];
    fname2   = [dir_where fname2suffix '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Aligning the sounds (looking at the maximum RMS value):
    opts.attenuate_by = attenuate_by;
    [signal1, signal2, fs] = il_wavread(fname1, fname2, opts);
    
    L = length(signal1);
    t = (1:L) / fs;

    % 3. Plotting time signals:
    fc = 20;
    yenv1 = il_get_envelope(signal1,fs,fc);
    yenv2 = il_get_envelope(signal2,fs,fc);

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

    % Copying the input files to the output folder:
    if bSave
        fname_out = sprintf('%s%s-att.wav',dir_out,fname1suffix);
        Wavwrite(signal1,fs,fname_out);
        fname_out = sprintf('%s%s-att.wav',dir_out,fname2suffix);
        Wavwrite(signal2,fs,fname_out);
    end
end

if bPart1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    method    = 3; % i = 3 = ERB;
    fname_out1 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname1suffix,method);
    fname_out2 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname2suffix,method);

    if bCreate
        [noise1 info4pede1 pede1] = icra5_noise4piano(signal1,fs,method);
        Wavwrite(noise1,fs,fname_out1);
        [noise2 info4pede2 pede2] = icra5_noise4piano(signal2,fs,method);
        Wavwrite(noise2,fs,fname_out2);
    else
        noise1  = Wavread(fname_out1);
        noise2  = Wavread(fname_out2);
    end

    N4fft = 1024; %8192*4;
    % signal1 = AM_random_noise(20,5000,77+3,2,44100);
        
    % Calculate LTAS
    [ltass1P, fHz,cal] = calcLTAS(signal1,fs,N4fft/44100);
    [ltass1N         ] = calcLTAS(noise1 ,fs,N4fft/44100);

    % Calculate LTAS
    [ltass2P, fHzP] = calcLTAS(signal2,fs,N4fft/44100);
    [ltass2N      ] = calcLTAS(noise2,fs,N4fft/44100);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    subplot(1,2,1);
    semilogx(   fHzP,ltass1P+cal,'-b', ... 
                fHzP,ltass1N+cal,'-k'); hold on, grid on;
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
    subplot(1,2,2);
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

    noise3 = 0.5*noise1+0.5*noise2;

    out3   = auditoryfilterbank(noise3,fs);
    RMSrel = rmsdb(out3)-max(rmsdb(out3));
    SNR    = 0;
    pede3  = icra_noise4piano_pedestal(noise3,fs,RMSrel,SNR);
    pede3 = pede3(1:length(noise3));
    fname_out3     = sprintf('%snoise-mixed-ICRA-meth-%.0f.wav'                 ,dir_out,method);
    fname_out3pede_p10_dB = sprintf('%snoise-mixed-ICRA-meth-%.0f+pede-SNR-%.0f-dB.wav',dir_out,method,SNR+10);
    fname_out3pede_0_dB   = sprintf('%snoise-mixed-ICRA-meth-%.0f+pede-SNR-%.0f-dB.wav',dir_out,method,SNR);
    fname_out3pede_m10_dB = sprintf('%snoise-mixed-ICRA-meth-%.0f+pede-SNR-%.0f-dB.wav',dir_out,method,SNR-10);

    noise_p10_dB = noise3+From_dB(-10)*pede3;
    noise_0_dB   = noise3+pede3;
    noise_m10_dB = noise3+From_dB(10)*pede3;

    if bSave
        Wavwrite(noise3      ,fs,fname_out3);
        Wavwrite(noise_p10_dB,fs,fname_out3pede_p10_dB);
        Wavwrite(noise_0_dB  ,fs,fname_out3pede_0_dB);
        Wavwrite(noise_m10_dB,fs,fname_out3pede_m10_dB);
    end
    
    figure;
    Text_ylabel = 'Amplitude';
    Text_xlabel = 'Time [s]';
    subplot(3,1,1)
    plot(t,noise1); hold on, grid on
    % plot(t,  nenv1,'k','LineWidth',2);
    % plot(t, -nenv1,'k','LineWidth',2);
    % plot(t_se1, From_dB(RMS_se1,10),'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('1');
    ha(end+1) = gca;

    subplot(3,1,2)
    plot(t, noise2,'r'); hold on, grid on
    % plot(t,  nenv2,'k','LineWidth',2);
    % plot(t, -nenv2,'k','LineWidth',2);
    % plot(t_se2, From_dB(RMS_se2,10),'k','LineWidth',2);
    ylabel(Text_ylabel)
    legend('2');
    ha(end+1) = gca;

    linkaxes(ha,'xy');
    xlabel(Text_xlabel)
    ylim([-0.049 0.049])
    xlim([0 max(t)])
    hFig(end+1) = gcf;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if bSave
        sig1_SNRm5 = signal1 + From_dB( 5) * noise_p10_dB;
        sig1_SNRm0 = signal1 + From_dB( 0) * noise_p10_dB;
        sig1_SNRp5 = signal1 + From_dB(-5) * noise_p10_dB;

        fname_out1 = sprintf('%s%s+mixed-noise-SNR-m05-dB-pede-p10dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRm5,fs,fname_out1);
        fname_out1 = sprintf('%s%s+mixed-noise-SNR-m00-dB-pede-p10dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRm0,fs,fname_out1);
        fname_out1 = sprintf('%s%s+mixed-noise-SNR-p05-dB-pede-p10dB.wav',dir_out,fname1suffix);
        Wavwrite(sig1_SNRp5,fs,fname_out1);
    end

    if bSave
        for i = 1:length(hFig)
            fname = sprintf('%sICRA-fig-%.0f',dir_out,i);
            Save_figure_as(hFig(i),fname,'epsc');
            Save_figure_as(hFig(i),fname,'fig');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    method    = 3; % i = 3 = ERB;
    out1_1 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname1suffix,method);
    out2_1 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname2suffix,method);
    [noise1_1 info4pede1 pede1_1] = icra5_noise4piano(signal1,fs,method);
    [noise1_2 info4pede1 pede1_2] = icra5_noise4piano(signal1,fs,method);
    [noise1_3 info4pede1 pede1_3] = icra5_noise4piano(signal1,fs,method);
    
    figure;
    plot(noise1_1)
    ha = gca;
    
    figure;
    plot(noise1_2)
    ha(end+1) = gca;
    
    figure;
    plot(noise1_3)
    ha(end+1) = gca;
    
    linkaxes(ha,'xy');
    
    fcl = 20;
    nenv1 = il_get_envelope(noise1_1,fs,fcl);
    nenv2 = il_get_envelope(noise1_2,fs,fcl);
    nenv3 = il_get_envelope(noise1_3,fs,fcl);
    
    figure;
    plot(nenv1); hold on
    plot(nenv2,'r--'); 
    plot(nenv3,'k-.'); 
    
end
   
if bPart3
    flow = 80;
    fhigh = 8000;
    fc = erbspacebw(flow, fhigh);
    
end

N = 8192*8;
K = N/2;
fs = 44100;
insig = [zeros(N/2-1,1); 1; zeros(N/2,1)]; % dirac

if bPart3
    
    [h_time,fc] = auditoryfilterbank(insig,fs);
    [h HdB f]   = freqfft2(h_time,K,fs);
    
    freqscut = [];
    fc_table = round(fc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get filterbank parameters:
    for i = 1:length(fc);

        f1 = Get_cutoff_freq_below_fc(f,abs(h(:,i)),fc(i));
        f2 = Get_cutoff_freq_above_fc(f,abs(h(:,i)),fc(i));
        f1_table = round(f1);
        f2_table = round(f2);
        freqscut = [freqscut; round(freqtoaud(fc(i))*10)/10 fc_table(i) f1_table f2_table];
        
    end
    BW = freqscut(:,4)-freqscut(:,3); % Bandwidth
    % BW(1) = freqscut(1,3); % To avoid not-a-number
    % estimation of centre frequency
    Q = transpose(fc)./BW;
    
    freqscut = [freqscut round(BW)]; 
    freqscut = [freqscut Round(Q,2)]; 
    
    var2latex(freqscut(1:16,:));
    var2latex(freqscut(17:31,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % figure;
    % plot(f,To_dB(abs(h))); grid on;
    % xlim([0 fmax2plot])
    % ylim([-18 5])
    % xlabel('Modulation frequency [Hz]')
    % ylabel('Attenuation [dB]')
    % title('Only modulation filterbank')
    
end

if bPart4
    
    opts = [];
    
    piano_1 = 'JBS50';
    take_1 = 3;
    
    piano_2 = 'NS19';
    take_2 = 1;
    
    attenuate_by = 10; % Because these signals were measured very closed to the soundboard (they are very loud)
    
    % 1. Sounds to be processed:
    fname1suffix = [piano_1 '-' note_test '_' num2str(take_1)];
    fname2suffix = [piano_2 '-' note_test '_' num2str(take_2)];
    fname1   = [dir_where fname1suffix '.wav'];
    fname2   = [dir_where fname2suffix '.wav'];
    
    % 2. Aligning the sounds (looking at the maximum RMS value):
    opts.attenuate_by = attenuate_by;
    [signal1, signal2, fs, s1, s2] = il_wavread(fname1, fname2, opts);
    [noise1 info4pede1 pede1] = icra_noise4piano(signal1,fs);
    [noise2 info4pede2 pede2] = icra_noise4piano(signal2,fs);

    % 2.1. Saving ICRA noise:
    fname_out1 = sprintf('%snoise-%s.wav',dir_out,fname1suffix);
    fname_out2 = sprintf('%snoise-%s.wav',dir_out,fname2suffix);
    Wavwrite(noise1,fs,fname_out1);
    Wavwrite(noise2,fs,fname_out2);
        
    L = length(signal1);
    t = (1:L)/fs; 
    
    %%%
    Lmax = max(length(s1),length(s2));
    Lmin = min(length(s1),length(s2));
    t1 = (1:length(s1))/fs;
    t2 = (1:length(s2))/fs;
    
    s1 = Do_cos_ramp(s1,fs,0,200);
    s2 = Do_cos_ramp(s2,fs,0,200);
    
    dur1 = length(s1)/fs;
    dur2 = length(s2)/fs;
    
    xoff1 = 0.175;
    xoff2 = 1.31;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N4fft = 1024; % for calculating the average spectrum:
    
    % s2 = Apply_IIR_Butter(s2,fs,50,'high');
    % s2 = Apply_IIR_Butter(s2,fs,8000,'low');
    [fHz, ltass1P, tRMS_s1, RMS_s1, cfModHzP, mtf1P] = il_get_signal_stats(    s1,fs,N4fft);
    [fHz, ltass1N, tRMS_s1, RMS_n1, cfModHzP, mtf1N] = il_get_signal_stats(noise1,fs,N4fft);
    [fHz, ltass2P, tRMS_s2, RMS_s2, cfModHzP, mtf2P] = il_get_signal_stats(    s2,fs,N4fft);
    [fHz, ltass2N, tRMS_s2, RMS_n2, cfModHzP, mtf2N] = il_get_signal_stats(noise2,fs,N4fft);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    leg1 = sprintf('Piano %s (take %.0f)',piano_1,take_1);
    leg2 = sprintf('Piano %s (take %.0f)',piano_2,take_2);
    
    % Figure 1:
    figure;
    plot(t2     , 2*s2); hold on; grid on
        
    legend(leg2)
    xlabel('Time [s]')
    ylabel('Amplitude [Pa]')
    hFig(1) = gcf;
    fname = [];
    fname{end+1} = [dir_out_figs sprintf('time-signal-%s-take-%.0f',piano_2,take_2)];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    % Figure 2:
    figure;
    subplot(1,2,1)
    plot(t2     , To_dB(2*s2/2e-5)); hold on; grid on
    plot(tRMS_s2,RMS_s2+100,'k','LineWidth',2);
    ha = gca;
    xlabel('Time [s]')
    ylabel('Sound Pressure Level [dB]')
    title(sprintf('Piano: %s-take-%.0f',piano_2,take_2))
    
    subplot(1,2,2)
    plot(t2     , To_dB(2*noise2/2e-5)); hold on; grid on
    plot(tRMS_s2,RMS_n2+100,'k','LineWidth',2);
    ha(end+1) = gca;
    % legend(leg1) 
    xlabel('Time [s]')
    ylabel('Sound Pressure Level [dB]')
    title('ICRA noise')
    hFig(end+1) = gcf;
    linkaxes(ha,'xy')
    ylim([45 75])
    xlim([0 dur2])
    fname{end+1} = [dir_out_figs sprintf('time-signal+noise-%s-take-%.0f',piano_2,take_2)];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    % Figure 3:
    figure;
    semilogx(   fHz,ltass1P,'-r'); hold on, grid on;
    semilogx(   fHz,ltass1N,'-k');
    xlabel('Frequency [Hz]')
    ylabel('Log. spectrum [dB]')
    legend({leg1,'noise (icra)'},'location','southwest')
    axis tight;
    XTick = [50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    hFig(end+1) = gcf;
    fname{end+1} = [dir_out_figs sprintf('LTAS-signal+noise-%s-take-%.0f',piano_1,take_1)];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    figure;
    semilogx(   fHz,ltass2P,'-b'); hold on, grid on;
    semilogx(   fHz,ltass2N,'-k');
    xlabel('Frequency [Hz]')
    ylabel('Log. spectrum [dB]')
    legend({leg2,'noise (icra)'},'location','southwest')
    axis tight;
    XTick = [50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTick);
    hFig(end+1) = gcf;
    fname{end+1} = [dir_out_figs sprintf('LTAS-signal+noise-%s-take-%.0f',piano_2,take_2)];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    % Figure 4:
    figure;
    plot(mtf2P,'.-','marker','s','markerSize',10,'color','b'), hold on
    plot(mtf2N,'.-','marker','x','markerSize',10,'color','k')
    set(gca,'xtick',1:numel(cfModHzP),'xticklabel',num2str(cfModHzP'))
    grid on;
    xlabel('Modulation frequency [Hz]')
    ylabel('MTF')
    legend({leg2,'noise (icra)'},'location','northwest')
    hFig(end+1) = gcf;
    fname{end+1} = [dir_out_figs sprintf('MTF-signal+noise-%s-take-%.0f',piano_2,take_2)];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    disp('')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figures 5:
    xa1 = ([0  dur2/dur1]+xoff1)/xoff2;
    xa2 = ([0          1]+xoff1)/xoff2;
    % Only for first subfigure:
    yoffset = 0.6;
    ya1 =  [0          0]+yoffset;
    ya2 =  [0.03    0.03]+yoffset;
     
    % Only for second subfigure:
    % yoffset = 0.15;
    % ya1_2 =  [0          0]+yoffset;
    % ya2_2 =  [0.03    0.03]+yoffset;
     
    %%%
    % Figure 5
    figure;
    plot(t1        , 2*s1,'r'); hold on
    plot(t1(1:Lmin), 2*s2); grid on

    xlim([0 dur1])

    annotation('doublearrow',xa1,ya1,'LineWidth',2);
    annotation('doublearrow',xa2,ya2,'LineStyle','--');
     
    legend(leg1,leg2)
    xlabel('Time [s]')
    ylabel('Amplitude [Pa]')
    hFig(end+1) = gcf;
    
    fname{end+1} = [dir_out_figs 'time-signals-short-long'];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    noise3 = 0.5*noise1 + 0.5*noise2;
        
    %%%
    figure
    subplot(1,3,1)
    plot(t,noise1,'r');
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    ha = gca;
    
    subplot(1,3,2)
    plot(t,noise2,'b');
    xlabel('Time [s]')
    ha(end+1) = gca;
    
    subplot(1,3,3)
    plot(t,noise3,'m'), hold on;
        
    xlabel('Time [s]')
    ha(end+1) = gca;
       
    linkaxes(ha,'xy');
    xlim([0 max(t)])
    ylim(0.045*[-1 1])
    
    opts.I_Width = 14;
    hFig(end+1) = Figure2paperfigureT(gcf,1,3,opts);
    
    fname{end+1} = [dir_out_figs 'time-ICRA-noises'];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    [RMS_n3 tRMS3] = rmsdb_sec(2*noise3,fs,10e-3,0.5);
    tRMS3 = tRMS3 + 10e-3/2;
    RMS_n3 = RMS_n3 + 100;
    
    [dB_max idx_max] = max(RMS_n3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SNR4pede = 10;
    out_tmp   = auditoryfilterbank(noise3,fs);
    RMSrel = rmsdb(out_tmp)-max(rmsdb(out_tmp));
    [pede noise3pede] = icra_noise4piano_pedestal(noise3,fs,RMSrel,SNR4pede);
    RMS_pede = rmsdb_sec(2*noise3pede,fs,10e-3,0.5);
    RMS_pede = RMS_pede + 100;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    subplot(1,3,1)
    plot(t,20*log10(abs(2*noise3)/2e-5),'m'); hold on
    plot(tRMS3, RMS_n3,'k','LineWidth',2);
    ylabel('Sound Pressure Level [dB]')
    xlabel('Time [s]')
    ha = gca;
    
    ofx = 0; %2;
    lline = 0.3; % length line in seconds
    subplot(1,3,2)
    plot(tRMS3, RMS_n3,'m','LineWidth',2); hold on
    plot([tRMS3(idx_max) tRMS3(idx_max)+lline],[dB_max dB_max],'k--','LineWidth',2)
    text(tRMS3(idx_max)+lline/3,dB_max+2,['L_m_a_x=' sprintf('%.1f dB',dB_max)],'FontSize',14)
    plot([0 max(t)],[dB_max-SNR4pede dB_max-SNR4pede],'k-.','LineWidth',2)
    text(   tRMS3(idx_max)+lline+ofx, dB_max-SNR4pede+0.5,['L_{pedestal}' sprintf('\n(SNR = %.0f dB)',SNR4pede)],'FontSize',14)
    
    xlabel('Time [s]')
    ha(end+1) = gca;
    
    ofx = 0; %-1;
    subplot(1,3,3)
    plot(tRMS3, RMS_pede,'LineWidth',2,'Color',rgb('Indigo')); hold on
    plot([tRMS3(idx_max) tRMS3(idx_max)+lline],[dB_max dB_max],'k--','LineWidth',2)
    text(tRMS3(idx_max)+lline/3,dB_max+2+ofx,['L_m_a_x'],'FontSize',14)
    plot([0 max(t)],[dB_max-SNR4pede dB_max-SNR4pede],'k-.','LineWidth',2);
    
    ofx = 0; %1;
    ofy = 0; %3;
    text(tRMS3(idx_max)+lline+ofx,dB_max-SNR4pede-8+ofy,['L_{pedestal}' sprintf('\n(SNR = %.0f dB)',SNR4pede)],'FontSize',14)
    
    xlabel('Time [s]')
    ha(end+1) = gca;
    set(ha,'FontSize',14);    
    linkaxes(ha,'xy')
    xlim([0 max(t)])
    ylim([35 75])
    
    hFig(end+1) = gcf; % Figure2paperfigureT(gcf,1,3,opts);
    
    fname{end+1} = [dir_out_figs 'time-ICRA-schematic'];
    if bSave
        Saveas(hFig(end),fname{end});
    end
    
    % get(gcf,'children')
    
    disp('')
    
end

rmpath([dir delim 'Tools'])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(['EOF: ' mfilename '.m'])

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal1, signal2, fs, sigorig1, sigorig2] = il_wavread(fname1,fname2,opts)
% function [signal1, signal2, fs, sigorig1, sigorig2] = il_wavread(fname1,fname2,opts)

if nargin < 3
    opts = [];
end
opts = Ensure_field(opts,'attenuate_by',0);
opts = Ensure_field(opts,'dur_ramp_down',150); % in ms, to be applied to truncated signal (the longer one)
opts = Ensure_field(opts,'window_length4max',10e-3); % in s
opts = Ensure_field(opts,'window_length4cal',0.5); % in s

attenuate_by = opts.attenuate_by;
dur_ramp_down = opts.dur_ramp_down;
window_length4max = opts.window_length4max;
window_length4cal = opts.window_length4cal;
% timeRMS2cal_bef = 150e-3; % assumes that the onset occurs at approx. t = 0.150 s
% timeRMS2cal_aft = 350e-3; % assumes that the target sound occurs in the next 0.350 s

[signal1, fs] = Wavread(fname1); 
[signal2, fs] = Wavread(fname2); 

t1 = ( 1:length(signal1) )/fs;
t2 = ( 1:length(signal2) )/fs;

[RMS_se1 t_se1] = rmsdb_sec(signal1,fs,window_length4max,0);
[RMS_se2 t_se2] = rmsdb_sec(signal2,fs,window_length4max,0);

%   2.1. Detecting the maximum and matching the levels:
[max_1 idx_1] = max(RMS_se1);
[max_2 idx_2] = max(RMS_se2);

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
% signal1 = Do_truncate(signal1,L);
% signal2 = Do_truncate(signal2,L);

RMS_se1 = rmsdb_sec(signal1,fs,window_length4cal);
RMS_se2 = rmsdb_sec(signal2,fs,window_length4cal);

max_1 = max(RMS_se1);
max_2 = max(RMS_se2);

delta_dB = max_2 - max_1;

if delta_dB >= 0 & delta_dB < 10 % then signal2 is louder (by less than 10 dB
    signal2 = From_dB(-abs(delta_dB)) * signal2; % attenuation of the louder signal
elseif abs(delta_dB) < 10
    signal1 = From_dB(-abs(delta_dB)) * signal1;
else
    error('The two signals you are trying to compare have very different levels')
end

signal1 = From_dB(-attenuate_by) * signal1;
signal2 = From_dB(-attenuate_by) * signal2;

if nargout > 3
    sigorig1 = signal1;
    sigorig2 = signal2;
end

signal1 = Do_truncate(signal1,L);
signal2 = Do_truncate(signal2,L);

if idxL == 1
    signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp_down); % then signal2 was truncated
else
    signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp_down); % then signal1 was truncated
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fHz, ltass, tRMS, RMS_dB, cfModHz, mtf] = il_get_signal_stats(insig,fs,N4fft)

[ltass, fHz, cal] = calcLTAS(insig,fs,N4fft/fs);
ltass = ltass + cal;

% RMS calculated every 10 ms, 50% overlap
[RMS_dB tRMS] = rmsdb_sec(2*insig,fs,10e-3,0.5);
tRMS = tRMS + 10e-3/2;

% MTF:
% cfModHz = [0.5 1 4 8 16 32];
[mtf,cfModHz] = calcMTF(insig,fs);