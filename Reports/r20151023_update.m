function y = r20151023_update(x)
% function y = r20151023_update(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 23/10/2015
% Last update on: 23/10/2015 
% Last use on   : 23/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

close all

dirout = [Get_TUe_paths('lx_Text') 'lx2015-10-23-decision-CASP' delim 'MATLAB' delim];
Mkdir(dirout);
diroutfigs = [dirout 'Figures' delim];
Mkdir(diroutfigs);
diroutaudio = [dirout 'Audio' delim];
Mkdir(diroutaudio);

dBFS = 100;
bSave = 0;

% bDoIOfunctions      = 0;
do_outermiddleear   = 0;
% bDoExcitation       = 0;
% bAnalysisTemplates  = 1;

fs = 44100;
ft2plotERB = 3:33; % ERB
ft2plot = audtofreq(ft2plotERB,'erb');

if do_outermiddleear
    
    N = 8192;
    K = N/2;
    h = [zeros(N/2,1); 1; zeros(N/2-1,1)];
    hp_fir = headphonefilter(fs); % Getting the filter coefficients at fs
    x1 = filter(hp_fir,1,h);
    ydBoe = 20*log10(abs(freqz(x1,1,K)));
    
    windowtype = 'hanning';
    me_fir   = middleearfilter(fs); % do_middleear - Lopez-Poveda
    me_fir_j = middleearfilter(fs,'jepsenmiddleear'); % do_middleear - Jepsen2008 (but not documented)
    
    xtmp = filter(me_fir,1,h);
    ydBme = 20*log10(abs(freqz(xtmp,1,K)));
    
    xtmp = filter(me_fir_j,1,h);
    ydBme_j = 20*log10(abs(freqz(xtmp,1,K)));
    
    xt   = filter(me_fir  ,1,x1);
    xt_j = filter(me_fir_j,1,x1);
    
    [xx xx f] = freqfft2(xt,K,fs,windowtype,dBFS,1); 
    ydBt = 20*log10(abs(freqz(xt,1,K)));
    ydBt_j = 20*log10(abs(freqz(xt_j,1,K)));
    
    calfactor = max(ydBt);
    calfactor2 = max(ydBt_j);
    
    figure;
    subplot(3,1,1)
    semilogx(f,ydBoe          ,'r','LineWidth',2); grid on, hold on
    semilogx(f,ydBme-calfactor,'b'); 
    semilogx(f,ydBt-+calfactor,'k--','LineWidth',2);
    legend('oe', ...
            sprintf('me - Lopez-Poveda + %.1f dB', abs(calfactor)), ...
            sprintf('oe + me + %.1f dB', abs(calfactor)));
    
    set(gca,'XTick',ft2plot(1:3:end));
    set(gca,'XTickLabel',round(ft2plot(1:3:end)));
    xlim([min(ft2plot) max(ft2plot)])
    
    ylabel('Gain [dB]')
    ylim([-35 28])
    xlabel('Frequency [Hz]')
    
    subplot(3,1,2)
    semilogx(f,ydBoe             ,'r','LineWidth',2); grid on, hold on
    semilogx(f,ydBme_j-calfactor2,'b'); 
    semilogx(f,ydBt_j-+calfactor2,'k--','LineWidth',2);
    legend('oe', ...
            sprintf('me - jepsen + %.1f dB', abs(calfactor2)), ...
            sprintf('oe + me + %.1f dB', abs(calfactor2)));
    
    set(gca,'XTick',ft2plot(1:3:end));
    set(gca,'XTickLabel',round(ft2plot(1:3:end)));
    xlim([min(ft2plot) max(ft2plot)])
    
    subplot(3,1,3)
    diffresps = (ydBt_j-+calfactor2) - (ydBt-+calfactor);
    semilogx(f,diffresps,'k--','LineWidth',2); grid on
    set(gca,'XTick',ft2plot(1:3:end));
    set(gca,'XTickLabel',round(ft2plot(1:3:end)));
    xlim([min(ft2plot) max(ft2plot)])
        
    ylabel('Difference [dB]')
    xlabel('Frequency [Hz]')
    legend('Jepsen - Lopez-Poveda')
    hFig1(end+1) = gcf;
    
    if bSave
        Saveas(hFig1(end),[diroutfigs 'om-ear']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bDoTFMF     = 0;
bDoTFDRNL   = 0;
bDoTFDRNL2  = 1; 
type = 'middleear'; % type = 'jepsenmiddleear';

h = [];
dire = [Get_TUe_paths('outputs') 'AMTControl-examples' delim]; % 'D:\Output\AMTControl-examples\'
dBFS = 100;

fmax2plot = 1000;
N = 8192*8;
K = N/2;
fs = 44100;
insig = [zeros(N/2-1,1); 1; zeros(N/2,1)]; % dirac

if bDoTFMF
    
    % To do: compare outsig with TFs
    fc_fig2b = 4000; % to get all the modulation filterbanks
    [outsig,mfc,params] = modfilterbank_debug(insig,fs,fc_fig2b);
    outsig = outsig{1};
    
    [x xdB f] = freqfft2(insig,K,fs);
    
    h(:,1) = freqz(params.b_mf1, params.a_mf1, K);
        
    for i = 2:length(mfc)
        
        exp1 = sprintf('h(:,%.0f) = freqz(params.b_mf%.0f, params.a_mf%.0f, K);',i,i,i);
        eval(exp1);
        
    end
    
    HdB = To_dB(abs(h));
    freqscut = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get filterbank parameters:
    for i = 1:length(mfc)
        
        if i ~= 1 
            f1 = Get_cutoff_freq_below_fc(f,abs(h(:,i)),mfc(i));
            f2 = Get_cutoff_freq_above_fc(f,abs(h(:,i)),mfc(i));
        else
            f1 = 0;
            f2 = Get_cutoff_freq_above_fc(f,abs(h(:,i)),mfc(i));
            mfc(1) = f2/2;
        end
        freqscut = [freqscut; mfc(i) f1 f2];
    end
    BW = freqscut(:,3)-freqscut(:,2); % Bandwidth
    BW(1) = freqscut(1,3); % To avoid not-a-number
     % estimation of centre frequency
    Q = transpose(mfc)./BW;
    
    freqscut = [freqscut BW]; 
    freqscut = [freqscut  Q]; % BW
    freqscut = Round(freqscut,2);
    var2latex(freqscut);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    plot(f,To_dB(abs(h))); grid on;
    xlim([0 fmax2plot])
    ylim([-18 5])
    xlabel('Modulation frequency [Hz]')
    ylabel('Attenuation [dB]')
    title('Only modulation filterbank')
    
    figure;
    hhigh = freqz(params.b_highest,params.a_highest,K);
    hhigh = repmat(hhigh,1,length(mfc));
    plot(f,To_dB(abs(h.*hhigh))); grid on; hold on
    xlim([0 fmax2plot])
    ylim([-18 5])
    xlabel('Modulation frequency [Hz]')
    ylabel('Attenuation [dB]')
    plot(f,To_dB(hhigh(:,1)),'k--','LineWidth',2)
     
    bands2plot = 12;
    for i = 1:bands2plot
        hout(:,i) = freqz(outsig(:,i),1,K);
    end
    
    hout_dB = To_dB(abs(hout(:,1:bands2plot)));
    figure;
%     plot(f,To_dB(abs(h(:,1:bands2plot))),'LineWidth',2); hold on
    plot(f,hout_dB(:,1:bands2plot)); hold on
    grid on;
    xlim([0 100])
    ylim([-18 5])
end

if bDoTFDRNL
    
    SPL = [40 70];
    LineWidths = [1 2];
    
    bands2plot = 31;
    
    f = 250;
    dur = N/fs;
    
    insig = Create_sin(f,dur,fs);
    
    for k = 1:length(SPL) 
        
        insigM = setdbspl( insig,SPL(k) );
        [x xdB f] = freqfft2(insigM,K,fs);

        [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs);

        for i = 1:bands2plot
            hout_drnl(:,i) = freqz(outsigdrnl(:,i),1,K);
            hout_lin(:,i) = freqz(paramsouts.outsiglin(:,i),1,K);
            hout_nlin(:,i) = freqz(paramsouts.outsignlin(:,i),1,K);
        end

        hout_dB_drnl = To_dB(abs(hout_drnl(:,1:bands2plot)));
        
        Max_values      = max( abs(hout_drnl) );
        
        hout_dB_drnl_norm = To_dB(hout_drnl./repmat(Max_values,length(hout_dB_drnl),1));
        
        figure(100);
        plot(f,hout_dB_drnl(:,1:3:bands2plot),'LineWidth',LineWidths(k)); hold on
        
        figure(101);
        plot(f,hout_dB_drnl_norm(:,1:3:bands2plot),'LineWidth',LineWidths(k)); hold on
        
    end
    figure(100)
    grid on;
    xlim([0 10000])
    ylim([0 SPL(2)+10])
    title('Impulse')
    
    figure(101)
    grid on;
    xlim([0 10000])
    ylim([-(SPL(2)+10) 0])
    title('Impulse-normalised')
    
end

if bDoTFDRNL2
    
    bFig2ab = 1;
    bFig2c = 0;
    
    SPL = [0 20:10:90 100]; % warning('temporal level')
    SPL_fig2c = [30 60 90]; % dB
    
    fc_fig2a = [ 250  500 800 1000 4000];
    bandidx = round( freqtoaud(fc_fig2a) )-2; % [5 9 12 14 25]; % 250, 500, 800, 1000, 4000 Hz respectively
    
    fc_fig2b = [1000 2400 4000 8000];
    dur = N/fs;
    
    %% Generation of Figure 2.a
    if bFig2ab
        for j = 1:length(fc_fig2a)
            insig = Create_sin(fc_fig2a(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs,type);    % DRNL
                [outsiggamma fcgamma]          = auditoryfilterbank(insigM,fs); % Gamma-tone
                outgamma = outsiggamma(:,bandidx(j)); % band centred at 1 kHz
                
                out = outsigdrnl(:,bandidx(j));
                lvl2agamma(j,k) = rmsdb(outgamma);
                lvl2a(j,k) = rmsdb(out);
            end       

        end

        dB_corr = dBFS;
        % Figure 2.a
        figure;
        subplot(1,2,1)
        plot(   SPL,lvl2a(1,:)+dB_corr,'k.--','LineWidth',1 ); hold on
        plot(   SPL,lvl2a(2,:)+dB_corr,'go-.','LineWidth',2 )
        plot(   SPL,lvl2a(3,:)+dB_corr,'r--') % 1k
        plot(   SPL,lvl2a(4,:)+dB_corr,'r>-.','LineWidth',2) % 1k
        plot(   SPL,lvl2a(5,:)+dB_corr,'b-', 'LineWidth',2 ); % 4k
        legend('250 Hz','500 Hz','800 Hz','1 kHz','4 kHz')
        title(['I/O functions, different CFs. DRNL (' type ')'])
        xlabel('Input level [dB SPL]')
        ylabel('DRNL output [dB re 1 m/s]')
        xlim([0 100])
        ylim([0 100])
        grid on
        set(gca,'YTick',[0:10:100])
        set(gca,'YTickLabel',[0:10:100])
        set(gca,'XTick',[0:10:100])
        set(gca,'XTickLabel',[0:10:100])
        
        % Figure 2.b
        for j = 1:4
            insig = Create_sin(fc_fig2b(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs,type);

                out = outsigdrnl(:,bandidx(4)); % band centred at 4 kHz
                lvl2b(j,k) = rmsdb(out);
            end       

        end

        % Figure 2.b
        subplot(1,2,2)
        plot(SPL,lvl2b(1,:),'r--'); hold on % 1 kHz
        plot(SPL,lvl2b(2,:),'k.--','LineWidth',1) % 2.4 kHz
        plot(SPL,lvl2b(3,:),'b-', 'LineWidth',2) % 4 kHz
        plot(SPL,lvl2b(4,:),'go-.','LineWidth',2) % 8 kHz
        legend('1 kHz','2.4 kHz','4 kHz','8 kHz')
        title(['CF = 4 kHz, different stim channels, DRNL (' type ')'])
        xlabel('Input level [dB SPL]')
        ylabel('DRNL output [dB re 1 m/s]')
        xlim([0 100])
        ylim([-100 0])
        grid on
        set(gca,'YTick',[-100:10:0])
        set(gca,'YTickLabel',[-100:10:0])
        set(gca,'XTick',[0:10:100])
        set(gca,'XTickLabel',[0:10:100])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 2.a gamma
        figure;
        subplot(1,2,1)
        plot(   SPL,lvl2a(1,:),'k.--','LineWidth',1 ); hold on
        plot(   SPL,lvl2a(2,:),'go-.','LineWidth',2 )
        plot(   SPL,lvl2a(3,:),'r--') % 1k
        plot(   SPL,lvl2a(4,:),'b-', 'LineWidth',2 ); % 4k
        legend('250 Hz','500 Hz','1 kHz','4 kHz')
        title(['I/O functions, different CFs. DRNL: (' type ')'])
        xlabel('Input level [dB SPL]')
        ylabel('DRNL output [dB re 1 m/s]')
        xlim([0 100])
        ylim([-100 0])
        grid on
        
        subplot(1,2,2)
        plot(   SPL,lvl2agamma(1,:),'k.--','LineWidth',1 ); hold on
        plot(   SPL,lvl2agamma(2,:),'go-.','LineWidth',2 )
        plot(   SPL,lvl2agamma(3,:),'r--') % 1k
        plot(   SPL,lvl2agamma(4,:),'b-', 'LineWidth',2 ); % 4k
        legend('250 Hz','500 Hz','1 kHz','4 kHz')
        title(['I/O functions, different CFs. DRNL (' type ')'])
        xlabel('Input level [dB SPL]')
        ylabel('Gammatone output [dB]')
        xlim([0 100])
        ylim([-100 0])
        grid on
        h(end+1) = gcf;
    end
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
