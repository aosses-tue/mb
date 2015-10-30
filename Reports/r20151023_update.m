function r20151023_update
% function r20151023_update
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

dirout = [Get_TUe_paths('lx_Text') 'lx2015-10-30-decision-CASP' delim 'MATLAB' delim];
Mkdir(dirout);
diroutfigs = [dirout 'Figures' delim];
Mkdir(diroutfigs);
% diroutaudio = [dirout 'Audio' delim]; Mkdir(diroutaudio);

dBFS = 100;
bSave = 1;
hFigs = [];

do_outermiddleear     = 0;
do_jepsen2008_fig_2a  = 1; 
do_jepsen2008_fig_2b  = 1;

bDoTFMF     = 0;
bDoTFDRNL   = 0;

% type = 'middleear'; % type = 'jepsenmiddleear';
% type = 'lopezpoveda2001';
type = 'jepsen2008';

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
    hFigs(end+1) = gcf;
    
    if bSave
        Saveas(hFigs(end),[diroutfigs 'om-ear']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmax2plot = 1000;
N = 22050;
K = N/2;
fs = 44100;

% Parameters for Figures 2.a and 2.b:
fc_fig2a= [250  500 800 1000 4000];
bandidx = round( freqtoaud(fc_fig2a) )-2; % [5 9 12 14 25]; % 250, 500, 800, 1000, 4000 Hz respectively
bandidx_fig2b = round( freqtoaud(4000) )-2;

dur     = N/fs;
SPL     = [0 20:10:90 100]; 

if do_jepsen2008_fig_2a
    
    %% Generation of Figure 2.a
    
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

    if strcmp(type,'lopezpoveda2001')
        dB_corr = 0;
    elseif strcmp(type,'jepsen2008')
        dB_corr = dBFS;
    end

    % Figure 2.a
    figure;
    % subplot(1,2,1)
    plot(   SPL,lvl2a(1,:)+dB_corr,'k.--','LineWidth',1 ); hold on
    plot(   SPL,lvl2a(2,:)+dB_corr,'go-.','LineWidth',2 )
    plot(   SPL,lvl2a(3,:)+dB_corr,'r--') % 1k
    plot(   SPL,lvl2a(4,:)+dB_corr,'r>-.','LineWidth',2) % 1k
    plot(   SPL,lvl2a(5,:)+dB_corr,'b-', 'LineWidth',2 ); % 4k
    legend( sprintf('%.0f Hz',fc_fig2a(1)), ... 
            sprintf('%.0f Hz',fc_fig2a(2)), ...
            sprintf('%.0f Hz',fc_fig2a(3)), ...
            sprintf('%.0f Hz',fc_fig2a(4)), ...
            sprintf('%.0f Hz',fc_fig2a(5)),'Location','NorthWest')
    title(['I/O functions, different CFs. DRNL (' type ')'])
    xlabel('Input level [dB SPL]')
    if strcmp(type,'lopezpoveda2001')
        ylabel('DRNL output [dB re 1 m/s]')
    elseif strcmp(type,'jepsen2008')
        ylabel('DRNL output [dB SPL]')
    end
        
    xlim([0 100])
    ylim([-100 0]+dB_corr)
    grid on
    set(gca,'YTick',[-100:10:0]+dB_corr)
    set(gca,'YTickLabel',[-100:10:0]+dB_corr)
    set(gca,'XTick',[0:10:100])
    set(gca,'XTickLabel',[0:10:100])
      
    hFigs(end+1) = gcf;
    
    if bSave
        Saveas(hFigs(end),[diroutfigs 'IO-functions-on-freq-' type],'epsc');
    end
end

if do_jepsen2008_fig_2b
    
    if strcmp(type,'lopezpoveda2001')
        dB_corr = 0;
    elseif strcmp(type,'jepsen2008')
        dB_corr = dBFS;
    end
    
    fc_fig2b = [1000 2400 4000 8000];
    
    for j = 1:length(fc_fig2b)
        insig = Create_sin(fc_fig2b(j),dur,fs);
        for k = 1:length(SPL) 

            insigM = setdbspl( insig,SPL(k) );

            [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs,type);

            out = outsigdrnl(:,bandidx_fig2b); % band centred at 4 kHz
            lvl2b(j,k) = rmsdb(out);
        end       

    end

    % Figure 2.b
    figure;
    % subplot(1,2,2)
    plot(SPL,lvl2b(1,:)+dB_corr,'r--'); hold on % 1 kHz
    plot(SPL,lvl2b(2,:)+dB_corr,'k.--','LineWidth',1) % 2.4 kHz
    plot(SPL,lvl2b(3,:)+dB_corr,'b-', 'LineWidth',2) % 4 kHz
    plot(SPL,lvl2b(4,:)+dB_corr,'go-.','LineWidth',2) % 8 kHz
    legend( sprintf('%.1f kHz',fc_fig2b(1)/1000), ...
        sprintf('%.1f kHz',fc_fig2b(2)/1000), ...
        sprintf('%.1f kHz',fc_fig2b(3)/1000), ...
        sprintf('%.1f kHz',fc_fig2b(4)/1000),'Location','NorthWest')
    
    title(['CF = 4 kHz, different stim channels, DRNL (' type ')'])
    xlabel('Input level [dB SPL]')
    if strcmp(type,'lopezpoveda2001')
        ylabel('DRNL output [dB re 1 m/s]')
    elseif strcmp(type,'jepsen2008')
        ylabel('DRNL output [dB SPL]')
    end
    xlim([0 100])
    ylim([-100+dB_corr 0+dB_corr])
    grid on
    set(gca,'YTick',[-100:10:0]+dB_corr)
    set(gca,'YTickLabel',[-100:10:0]+dB_corr)
    set(gca,'XTick',[0:10:100])
    set(gca,'XTickLabel',[0:10:100])
    
    hFigs(end+1) = gcf;
    
    if bSave
        Saveas(hFigs(end),[diroutfigs 'IO-functions-off-freq-' type],'epsc');
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
