function y = r20151016_update_drnl
% function y = r20151016_update_drnl
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20150814_update
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original file : r20150814_update.m
% Created on    : 12/08/2015
% Last update on: 15/08/2015 
% Last use on   : 15/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

h = [];
dire = [Get_TUe_paths('outputs') 'AMTControl-examples' delim]; % 'D:\Output\AMTControl-examples\'
dBFS = 100;

fmax2plot = 1000;
N = 8192*8;
K = 22050;
fs = 44100;

bFig2ab = 1;
bFig2c = 0;

do_lopezpoveda_fig3 = 1;
do_lopezpoveda_fig5 = 0;
do_jepsen_fig2ab    = 0;
do_jepsen_fig2c     = 1;

if do_lopezpoveda_fig3
    
    bUse_drnl_from_AMT = 0;
    bUse_drnl_own = ~bUse_drnl_from_AMT;
    
    kv.predrnl = {};
    kv.postdrnl = {};
    
	f1000_debug = {'flow',1000, 'fhigh',1000};
    
    %% Lopez-Poveda and Meddis 2001, Figure 3, b)
    % input signal: 50ms pure tones, sampled at 44.1kHz
    fs = 44100;
    T = 0.05;       
    t = (0:1/fs:T - 1/fs)';
    fsig = 250:25:1750;    

    result3b  = zeros(1,length(fsig));
    lin3b     = zeros(1,length(fsig));
    nlin3b    = zeros(1,length(fsig));

    level = 20e-6 * 10^(30/20);
  
    for ii = 1:length(fsig)

        insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  

        if bUse_drnl_from_AMT
            hp_fir = headphonefilter(fs);
            insig = filter(hp_fir,1,insig);
            [ y_lin, ~]     = drnl(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:});    
            [y_nlin, ~]     = drnl(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:});
        end
        
        if bUse_drnl_own
            insig = gaindb(insig,-6); % own calibration
            [ y_lin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:},'no_output_gain');    
            [y_nlin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:},'no_output_gain');
        end
        outsig          = y_lin + y_nlin;

        result3b(1,ii)  = rms(outsig(floor(length(insig)/2):end));
        lin3b(1,ii)     = rms(y_lin(floor(length(insig)/2):end));
        nlin3b(1,ii)    = rms(y_nlin(floor(length(insig)/2):end));
        
    end
    %%% 3c
    result3c  = zeros(1,length(fsig));
    lin3c     = zeros(1,length(fsig));
    nlin3c    = zeros(1,length(fsig));

    level = 20e-6 * 10^(85/20);
  
    for ii = 1:length(fsig)

        insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  
        
        if bUse_drnl_from_AMT
            hp_fir = headphonefilter(fs);
            insig  = filter(hp_fir,1,insig);
            [ y_lin, ~]     = drnl(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:});    
            [y_nlin, ~]     = drnl(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:});
        end
        if bUse_drnl_own
            insig = gaindb(insig,-6); % own calibration
            [ y_lin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:},'no_output_gain');    
            [y_nlin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:},'no_output_gain');
        end
        outsig          = y_lin + y_nlin;

        result3c(1,ii)  = rms(outsig(floor(length(insig)/2):end));
        lin3c(1,ii)     = rms(y_lin(floor(length(insig)/2):end));
        nlin3c(1,ii)    = rms(y_nlin(floor(length(insig)/2):end));
        
    end
        
    figure
    set(gcf,'Position',[50,50,400,760])
    subplot(2,1,1)
    plot(fsig,result3b)
    hold on
    plot(fsig,lin3b,'-.g')
    plot(fsig,nlin3b,':r')
    set(gca,'YScale','log')
    % grid on
    set(gca,'XLim',[250 1750],'Layer','top')
    set(gca,'YLim',[1e-07 1e-03],'Layer','top')
    set(gca,'Position',[0.285,0.5838,0.62,0.3412]);
    title('30 dB SPL')
    xlabel('Frequency (Hz)')
    ylabel('DRNL filter output (m/s)')

    subplot(2,1,2)
    plot(fsig,result3c)
    hold on
    plot(fsig,lin3c,'-.g')
    plot(fsig,nlin3c,':r')
    set(gca,'YScale','log')
    % grid on
    set(gca,'XLim',[250 1750],'Layer','top')
    set(gca,'YLim',[1e-05 1e-01],'Layer','top')
    set(gca,'Position',[0.285,0.11,0.62,0.3412]);
    title('85 dB SPL')
    xlabel('Frequency (Hz)')
    ylabel('DRNL filter output (m/s)')
    leg=legend('DRNL output', 'Linear path output', 'Nonlinear path output');
    set(leg,'Position',[0.0133, 0.4759, 0.4525, 0.0798]);
  
end

if do_lopezpoveda_fig5
    
    ftest = [250 500 800 1000 1500 2000 4000 8000];
    
    definput.import        ={'drnl_CASP' };
    definput.importdefaults={'lopezpoveda2001'};
    definput.keyvals.subfs=[];
    definput.keyvals.subfs=[];

    [flags,keyvalsLP2001]  = ltfatarghelper({'flow','fhigh'},definput,{});

    oCFlin = il_polfun(keyvalsLP2001.lin_fc,ftest);
    oBWlin = il_polfun(keyvalsLP2001.lin_bw,ftest);
    oCFnl  = il_polfun(keyvalsLP2001.nlin_fc_before,ftest);
    oBWnl  = il_polfun(keyvalsLP2001.nlin_bw_before,ftest);
    oa  = il_polfun(keyvalsLP2001.nlin_a,ftest);
    ob  = il_polfun(keyvalsLP2001.nlin_b,ftest);
    oc  = il_polfun(keyvalsLP2001.nlin_c,ftest);
    og  = il_polfun(keyvalsLP2001.lin_gain,ftest);
        
    colors = {'b-','r-','g-','k-','bo-','ro-','go-','ko-'};
    legends = {'CFlin','BWlin','CFnl','BWnl','a','b','c','d'};
    
    figure;
    subplot(1,2,1)
    loglog(ftest,oCFlin, colors{1}); hold on; grid on
    plot(ftest,oBWlin,colors{2});
    plot(ftest,oCFnl,colors{3});
    plot(ftest,oBWnl,colors{4});
    plot(ftest,oa,colors{5});
    plot(ftest,ob,colors{6});
    plot(ftest,oc,colors{7});
    plot(ftest,og,colors{8});
    
    legend(legends,'Location','NorthWest');
    
    definput.importdefaults={'jepsen2008'};
    [flags,keyvalsJ2008]  = ltfatarghelper({'flow','fhigh'},definput,{});
    
    oCFlin = il_polfun(keyvalsJ2008.lin_fc,ftest);
    oBWlin = il_polfun(keyvalsJ2008.lin_bw,ftest);
    oCFnl  = il_polfun(keyvalsJ2008.nlin_fc_before,ftest);
    oBWnl  = il_polfun(keyvalsJ2008.nlin_bw_before,ftest);
    fabove = 1500;
    oa_above1500 = il_polfun(keyvalsJ2008.nlin_a_above,fabove*ones(size(ftest)));
    ob_above1500 = il_polfun(keyvalsJ2008.nlin_b_above,fabove*ones(size(ftest)));
    oc_above1500 = il_polfun(keyvalsJ2008.nlin_c,fabove*ones(size(ftest)));
    
    oa  = il_polfun(keyvalsJ2008.nlin_a,ftest);
    ob  = il_polfun(keyvalsJ2008.nlin_b,ftest);
    oc  = il_polfun(keyvalsJ2008.nlin_c,ftest);
    og  = il_polfun(keyvalsJ2008.lin_gain,ftest);
    
    idxabove = find(ftest >= 1500);
    oa(idxabove) = oa_above1500(idxabove);
    ob(idxabove) = ob_above1500(idxabove);
    oc(idxabove) = oc_above1500(idxabove);
    
    ha = gca;
    
    subplot(1,2,2)
    loglog(ftest,oCFlin, colors{1}); hold on; grid on
    plot(ftest,oBWlin,colors{2});
    plot(ftest,oCFnl,colors{3});
    plot(ftest,oBWnl,colors{4});
    plot(ftest,oa,colors{5});
    plot(ftest,ob,colors{6});
    plot(ftest,oc,colors{7});
    plot(ftest,og,colors{8});
    
    ha(end+1) = gca;
    linkaxes(ha,'xy')
    
end

if do_jepsen_fig2ab
    
    SPL = [0 20:10:90 100];
    
    % fc_fig2a = [100 250  500 1000 4000 8000];
    % fc_fig2a = [86.9 257.4 520 1055.8  3982.6 7819.2];
    fc_fig2a = [257.4 520 1055.8  3982.6];
    bandidxa = round(freqtoaud(fc_fig2a))-2;
    
    fc_fig2b = [1000 2400 4000 8000];
    bandidxb = round(freqtoaud(4000))-2;
    % fc_fig2c = Get_OB_freqs(3,250,2000);
    dur = N/fs;
    
    % mode = 'nlinonly';
    mode = 'bothparts';
    
    %% Generation of Figure 2.a
    if bFig2ab
        for j = 1:length(fc_fig2a)
            insig = Create_sin(fc_fig2a(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                % [outsig_lin fcs oparams] = drnl_CASP_debug(insigM,fs,'linonly');    % DRNL
                % [outsig_nlin fcs]        = drnl_CASP_debug(insigM,fs,'nlinonly');    % DRNL
                % out_lin     = outsig_lin(:,bandidx(j));
                % out_nlin    = outsig_nlin(:,bandidx(j));
                % out_tot = out_lin + out_nlin;

                [outsig fcs] = drnl_CASP_debug(insigM,fs,mode,'noouterear');    % DRNL
                out_tot = outsig(:,bandidxa(j));
                % lvl2agamma(j,k) = rmsdb(outgamma);

                Ns = round(length(insigM)/2);
                lvl2a(j,k)      = rmsdb(out_tot(Ns+1:end));
                
                disp('')
            end       

        end

        figure;
        % Figure 2.a
        subplot(1,2,1)
        plot(   SPL,lvl2a(1,:),'k.--','LineWidth',1 ); hold on
        plot(   SPL,lvl2a(2,:),'go-.','LineWidth',2 )
        plot(   SPL,lvl2a(3,:),'r--') % 1k
        plot(   SPL,lvl2a(4,:),'b-', 'LineWidth',2 ); % 4k
        % plot(   SPL,lvl2a(5,:),'m>-', 'LineWidth',2 ); % 4k
        
        legend( sprintf('%.0f Hz',fc_fig2a(1)), ... 
                sprintf('%.0f Hz',fc_fig2a(2)), ...
                sprintf('%.0f Hz',fc_fig2a(3)), ...
                sprintf('%.0f Hz',fc_fig2a(4)) ) 
        title('I/O functions, different CFs')
        xlabel('Input level [dB SPL]')
        ylabel('DRNL output [dB re 1 m/s]')
        xlim([0 100])
        ylim([-100 0])
        grid on
        set(gca,'YTick',[-100:10:0])
        set(gca,'YTickLabel',[-100:10:0])
        set(gca,'XTick',[0:10:100])
        set(gca,'XTickLabel',[0:10:100])
        
        % Figure 2.b
        for j = 1:length(fc_fig2b)
            insig = Create_sin(fc_fig2b(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                [outsig fcs oparams] = drnl_CASP_debug(insigM,fs,mode,'noouterear');

                out_tot = outsig(:,bandidxb); % band centred at 4 kHz
                Ns = round(length(insigM)/2);
                lvl2b(j,k) = rmsdb(out_tot(Ns+1:end));
                
            end       

        end

        % Figure 2.b
        subplot(1,2,2)
        plot(SPL,lvl2b(1,:),'r--'); hold on % 1 kHz
        plot(SPL,lvl2b(2,:),'k.--','LineWidth',1) % 2.4 kHz
        plot(SPL,lvl2b(3,:),'b-', 'LineWidth',2) % 4 kHz
        plot(SPL,lvl2b(4,:),'go-.','LineWidth',2) % 8 kHz
        
        legend( sprintf('%.1f kHz',fc_fig2b(1)/1000), ...
                sprintf('%.1f kHz',fc_fig2b(2)/1000), ...
                sprintf('%.1f kHz',fc_fig2b(3)/1000), ...
                sprintf('%.1f kHz',fc_fig2b(4)/1000))
        
        title('CF = 4 kHz, different stim channels')
        xlabel('Input level [dB SPL]')
        ylabel('DRNL output [dB re 1 m/s]')
        xlim([0 100])
        ylim([-100 0])
        grid on
        set(gca,'YTick',[-100:10:0])
        set(gca,'YTickLabel',[-100:10:0])
        set(gca,'XTick',[0:10:100])
        set(gca,'XTickLabel',[0:10:100])
        
        h(end+1) = gcf;
        
    end
    
end

if do_jepsen_fig2c
    
    kv.predrnl = {};
    kv.postdrnl = {};
    
	f1000_debug = {'flow',1000, 'fhigh',1000};
    
    fs = 44100;
    T = 0.05;       
    t = (0:1/fs:T - 1/fs)';
    fsig = 250:25:1750;    

    result2c  = zeros(1,length(fsig));
    lin2c     = zeros(1,length(fsig));
    nlin2c    = zeros(1,length(fsig));

    level = 20e-6 * 10^(30/20);
  
    for ii = 1:length(fsig)

        insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  

        insig = gaindb(insig,-6); % own calibration
        [ y_lin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:});    
        [y_nlin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:});
        
        outsig          = y_lin + y_nlin;

        result2c(1,ii)  = rmsdb(outsig(floor(length(insig)/2):end));
        lin2c(1,ii)     = rmsdb(y_lin(floor(length(insig)/2):end));
        nlin2c(1,ii)    = rmsdb(y_nlin(floor(length(insig)/2):end));
        
    end
    %%% 2e
    result2d  = zeros(1,length(fsig));
    lin2d     = zeros(1,length(fsig));
    nlin2d    = zeros(1,length(fsig));

    level = 20e-6 * 10^(60/20);
  
    for ii = 1:length(fsig)

        insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  
        
        insig = gaindb(insig,-6); % own calibration
        [ y_lin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:});    
        [y_nlin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:});
        outsig          = y_lin + y_nlin;

        result2d(1,ii)  = rmsdb(outsig(floor(length(insig)/2):end));
        lin2d(1,ii)     = rmsdb(y_lin(floor(length(insig)/2):end));
        nlin2d(1,ii)    = rmsdb(y_nlin(floor(length(insig)/2):end));
        
    end
     %%% 3c
    result2e  = zeros(1,length(fsig));
    lin2e     = zeros(1,length(fsig));
    nlin2e    = zeros(1,length(fsig));

    level = 20e-6 * 10^(60/20);
  
    for ii = 1:length(fsig)

        insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  
        
        insig = gaindb(insig,-6); % own calibration
        [ y_lin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'linonly', kv.postdrnl{:});    
        [y_nlin, ~]     = drnl_CASP_debug(insig, fs, kv.predrnl{:},f1000_debug{:},'nlinonly', kv.postdrnl{:});
        outsig          = y_lin + y_nlin;

        result2e(1,ii)  = rmsdb(outsig(floor(length(insig)/2):end));
        lin2e(1,ii)     = rmsdb(y_lin(floor(length(insig)/2):end));
        nlin2e(1,ii)    = rmsdb(y_nlin(floor(length(insig)/2):end));
        
    end   
    figure
    set(gcf,'Position',[50,50,400,760])
    subplot(3,1,1)
    plot(fsig,result2c)
    hold on
    plot(fsig,lin2c,'-.g')
    plot(fsig,nlin2c,':r')
    set(gca,'XLim',[250 1750],'Layer','top')
    set(gca,'YLim',[-100 0],'Layer','top')
    % set(gca,'Position',[0.285,0.5838,0.62,0.3412]);
    title('30 dB SPL')
    xlabel('Frequency (Hz)')
    ylabel('DRNL filter output [dB]')

    subplot(3,1,2)
    plot(fsig,result2d)
    hold on
    plot(fsig,lin2d,'-.g')
    plot(fsig,nlin2d,':r')
    % grid on
    set(gca,'XLim',[250 1750],'Layer','top')
    set(gca,'YLim',[-100 0],'Layer','top')
    % set(gca,'Position',[0.285,0.11,0.62,0.3412]);
    title('60 dB SPL')
    xlabel('Frequency (Hz)')
    ylabel('DRNL filter output [dB]')
    leg=legend('DRNL output', 'Linear path output', 'Nonlinear path output');
    set(leg,'Position',[0.0133, 0.4759, 0.4525, 0.0798]);
    
    subplot(3,1,3)
    plot(fsig,result2e)
    hold on
    plot(fsig,lin2e,'-.g')
    plot(fsig,nlin2e,':r')
    % grid on
    set(gca,'XLim',[250 1750],'Layer','top')
    set(gca,'YLim',[-100 0],'Layer','top')
    % set(gca,'Position',[0.285,0.11,0.62,0.3412]);
    title('90 dB SPL')
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(['EOF: ' mfilename '.m'])
% Inline functions:

function outpar=il_polfun(par,fc)
% p0 = par(1)
% m  = par(2)
% outpar=10^(par(1)+par(2)*log10(fc));
outpar=10^(par(1))*fc.^par(2);

disp('')
