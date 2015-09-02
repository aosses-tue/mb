function r20150901_update(typetone,SPL)
% function r20150901_update(typetone,SPL)
%
% 1. Description:
%       Computing excitation patterns of the DRNL
% 
% 2. Stand-alone example:
%       typetone = 1; % sine tones
%       SPL = 60; % dB SPL
%       r20150901_update(typetone,SPL)
% 
%       typetone = 2; % 20-Hz Gaussian noise
%       SPL = 60; % dB SPL
%       r20150901_update(typetone,SPL)
%
%       typetone = 3; % 100-Hz Gaussian noise
%       SPL = 60; % dB SPL
%       r20150901_update(typetone,SPL)
%
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 27/08/2015
% Last update on: 01/09/2015 
% Last use on   : 01/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    typetone = 1; % sinusoids
end

if nargin < 2
    SPL = 60;
end

bDiary = 0;
Diary(mfilename,bDiary);
close all

bSave = 1;
hFig1 = [];
hFig2 = [];
hFig3 = [];

dirout = [Get_TUe_paths('lx_Text') 'lx2015-09-01-update' delim 'MATLAB' delim];
Mkdir(dirout);
diroutfigs = [dirout 'Figures' delim];
Mkdir(diroutfigs);
diroutaudio = [dirout 'Audio' delim];
Mkdir(diroutaudio);


dBFS = 100;

bDoIOfunctions      = 0;
do_outermiddleear   = 0;
bDoExcitation       = 0;
bAnalysisTemplates  = 1;

fs = 44100;
ft2plotERB = 3:33; % ERB
ft2plot = audtofreq(ft2plotERB,'erb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoIOfunctions
    r20150814_update([0 0 0 0 1]); % then store figures manually
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_outermiddleear
    
    N = 8192;
    K = N/2;
    h = [zeros(N/2,1); 1; zeros(N/2-1,1)];
    hp_fir = headphonefilter(fs); % Getting the filter coefficients at fs
    x1 = filter(hp_fir,1,h);
    ydBoe = 20*log10(abs(freqz(x1,1,K)));
    
    windowtype = 'hanning';
    % [y ydBoe f] = freqfft2(x1,K,fs,windowtype,dBFS,1); 
    me_fir = middleearfilter(fs); % do_middleear - Lopez-Poveda
    
    x2 = filter(me_fir,1,h);
    % [y ydBme f] = freqfft2(x2,K,fs,windowtype,dBFS,1); 
    ydBme = 20*log10(abs(freqz(x2,1,K)));
    
    xt = filter(me_fir,1,x1);
    [xx xx f] = freqfft2(xt,K,fs,windowtype,dBFS,1); 
    ydBt = 20*log10(abs(freqz(xt,1,K)));
    
    calfactor = max(ydBt);
    
    figure;
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
    
    hFig1(end+1) = gcf;
    
    if bSave
        Saveas(hFig1(end),[diroutfigs 'om-ear']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExcitation
    
    fton    = 1300;
    ftoff   = 2000;
    [xx, idxon ] = min(abs(ft2plot - fton ));
    [xx, idxoff] = min(abs(ft2plot - ftoff));
    SPL = 60; % :6:84;
    dur = 1; % s

    nfc = length(fton);
    
    nSPL     = length(SPL);
    lvl_drnl = nan(nfc,length(ft2plotERB));
    lvltot   = nan(nfc,1);

    lvl_gt = nan(nfc,length(ft2plotERB));
    lvltot_gt = nan(nfc,1);
    step_bands = 1; % set to 1 to analyse every band

    insigS = [];
    insigBW100 = [];
    insigBW20 = [];
    ha1 = [];
    ha2 = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:step_bands:nfc

        switch typetone
            case 1 % sinusoid
                insig = Create_sin(fton(j),dur,fs);
                insigS = [insigS insig];
                tonelabel = 'Test tone';
                tonelabel_short = 'T';
            case 2 % 20-Hz Gaussian
                BW = 20;
                insig = AM_random_noise_BW(fton(j),BW,60,dur,fs);
                insigBW20 = [insigBW20 insig];
                tonelabel = '20-Hz GN signal';
                tonelabel_short = 'GN-020-Hz';
                W(j) = Get_envelope_descriptors(insig,'W');
                V(j) = Get_envelope_descriptors(insig,'V');
            case 3 % 100-Hz Gaussian
                BW = 100;
                insig = AM_random_noise_BW(fton(j),BW,60,dur,fs);
                insigBW100 = [insigBW100 insig];
                tonelabel = '100-Hz GN signal';
                tonelabel_short = 'GN-100-Hz';
                W(j) = Get_envelope_descriptors(insig,'W');
                V(j) = Get_envelope_descriptors(insig,'V');
        end

        for k = 1:nSPL 

            insigM = setdbspl( insig,SPL(k) );

            [out_drnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs); 
            outsig_wav_tot      = sum(out_drnl,2);
            outsig_wav_onfreq   = out_drnl(:,idxon );
            outsig_wav_offfreq  = out_drnl(:,idxoff);
            
            Wavwrite(outsig_wav_tot    ,fs,sprintf('%sdrnl-%s-fc-%.0f-Hz-%.0f-dB'                ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_onfreq ,fs,sprintf('%sdrnl-%s-fc-%.0f-Hz-%.0f-dB-onfreq'         ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_offfreq,fs,sprintf('%sdrnl-%s-fc-%.0f-Hz-%.0f-dB-offfreq-%.0f-Hz',diroutaudio,tonelabel_short,fton,SPL(k),ftoff))
            
            lvl_drnl(j,:) = rmsdb(out_drnl) + dBFS; % level per band
            out_drnl_t(:,j) = outsig_wav_tot;
            out_drnl_on_t(:,j) = outsig_wav_onfreq;
            out_drnl_off_t(:,j) = outsig_wav_offfreq;
            
            lvltotdrnl(j,1) = sum_db( rmsdb(out_drnl) ) + dBFS;
            lvltotdrnl(j,2) = sum_db( rmsdb(paramsouts.outsiglin) ) + dBFS;
            lvltotdrnl(j,3) = sum_db( rmsdb(paramsouts.outsignlin) ) + dBFS;
            lvltotdrnl(j,4) = 100*sum(rms(paramsouts.outsignlin)) / sum(rms(out_drnl));

            ingamma = gaindb(paramsouts.insig_after_TF,-paramsouts.me_gain_TF);
            % ingamma = insigM;
            [out_gt fcgamma]= auditoryfilterbank(ingamma,fs);
            outsig_wav_tot      = sum(out_gt,2);
            outsig_wav_onfreq   = out_gt(:,idxon );
            outsig_wav_offfreq  = out_gt(:,idxoff);
            
            Wavwrite(outsig_wav_tot    ,fs,sprintf('%sgt-%s-fc-%.0f-Hz-%.0f-dB'                ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_onfreq ,fs,sprintf('%sgt-%s-fc-%.0f-Hz-%.0f-dB-onfreq'         ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_offfreq,fs,sprintf('%sgt-%s-fc-%.0f-Hz-%.0f-dB-offfreq-%.0f-Hz',diroutaudio,tonelabel_short,fton,SPL(k),ftoff))
            
            lvl_gt(j,:)     = rmsdb(out_gt) + dBFS;
            lvltot_gt(j,1)  = sum_db( rmsdb(out_gt) ) + dBFS;

            out_gt_t(:,j) = outsig_wav_tot;
            out_gt_on_t(:,j)  = outsig_wav_onfreq;
            out_gt_off_t(:,j) = outsig_wav_offfreq;

            figure; 
            plot(1:31, lvl_drnl(j,:),'LineWidth',2), grid on, hold on
            plot(1:31, lvl_gt(j,:),'r')
            title(sprintf('Excitation pattern in DRNL filterbank (outer- and middle - ear compensation)\n%s centred at fc = %.1f [Hz]',tonelabel,fton(j)));
            legend('DRNL','4th-order gt','Location','SouthEast')

            xlabel('Frequency [Hz]')
            ylabel('Log-spectrum at DRNL output [dB]')

            set(gca,'XTick',[1:4:31]);
            set(gca,'XTickLabel',round(audtofreq([1:4:31]+2)));

            ha1(end+1) = gca;
            hFig2(end+1) = gcf;
        end       

    end

    linkaxes([ha1],'xy');

    switch typetone
        case 1
            ylim([SPL-95 SPL+5])
        case 2
            ylim([SPL-95 SPL+15])
        otherwise
            ylim([SPL-95 SPL+5])
    end

    xlim([0.5 31+0.5])

    if bSave
        for i = 1:length(hFig2) 
            Saveas(hFig2(i),sprintf('%sfig-DRNL-with-om-TF-%.0f-SPL-%.0f-dB-%s',diroutfigs,i,SPL,tonelabel_short)); 
        end
        % close all
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:step_bands:nfc

        switch typetone
            case 1 % sinusoid
                insig = Create_sin(fton(j),dur,fs);
            case 2 % 20-Hz Gaussian
                BW = 20;
                insig = insigBW20(:,j);
            case 3 % 100-Hz Gaussian
                BW = 100;
                insig = insigBW100(:,j);
        end

        for k = 1:nSPL 

            insigM = setdbspl( insig,SPL(k) );

            [out_drnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs,'nomiddleear','noouterear'); 
            outsig_wav_tot      = sum(out_drnl,2);
            outsig_wav_onfreq   = out_drnl(:,idxon );
            outsig_wav_offfreq  = out_drnl(:,idxoff);
            
            Wavwrite(outsig_wav_tot    ,fs,sprintf('%sdrnlnoear-%s-fc-%.0f-Hz-%.0f-dB'                ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_onfreq ,fs,sprintf('%sdrnlnoear-%s-fc-%.0f-Hz-%.0f-dB-onfreq'         ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_offfreq,fs,sprintf('%sdrnlnoear-%s-fc-%.0f-Hz-%.0f-dB-offfreq-%.0f-Hz',diroutaudio,tonelabel_short,fton,SPL(k),ftoff))
            
            lvl_drnl(j,:) = rmsdb(out_drnl) + dBFS; % level per band
            out_drnl_noear_t(:,j) = outsig_wav_tot;
            out_drnl_noear_on_t(:,j) = outsig_wav_onfreq;
            out_drnl_noear_off_t(:,j) = outsig_wav_offfreq;

            lvldrnl_noear(j,1) = sum_db( rmsdb(out_drnl) ) + dBFS;
            lvldrnl_noear(j,2) = sum_db( rmsdb(paramsouts.outsiglin) ) + dBFS;
            lvldrnl_noear(j,3) = sum_db( rmsdb(paramsouts.outsignlin) ) + dBFS;
            lvldrnl_noear(j,4) = 100*sum(rms(paramsouts.outsignlin)) / sum(rms(out_drnl));

            
            ingamma = gaindb(paramsouts.insig_after_TF,-paramsouts.me_gain_TF);
            % ingamma = insigM;
            [out_gt fcgamma]= auditoryfilterbank(ingamma,fs);
            
            outsig_wav_tot      = sum(out_gt,2);
            outsig_wav_onfreq   = out_gt(:,idxon );
            outsig_wav_offfreq  = out_gt(:,idxoff);
            
            Wavwrite(outsig_wav_tot    ,fs,sprintf('%sgtnoear-%s-fc-%.0f-Hz-%.0f-dB'                ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_onfreq ,fs,sprintf('%sgtnoear-%s-fc-%.0f-Hz-%.0f-dB-onfreq'         ,diroutaudio,tonelabel_short,fton,SPL(k)))
            Wavwrite(outsig_wav_offfreq,fs,sprintf('%sgtnoear-%s-fc-%.0f-Hz-%.0f-dB-offfreq-%.0f-Hz',diroutaudio,tonelabel_short,fton,SPL(k),ftoff))
            
            lvl_gt(j,:)     = rmsdb(out_gt) + dBFS;
            lvltot_gt(j,1)  = sum_db( rmsdb(out_gt) ) + dBFS;

            out_gt_noear_t(:,j) = outsig_wav_tot;
            out_gt_noear_on_t(:,j) = outsig_wav_onfreq;
            out_gt_noear_off_t(:,j) = outsig_wav_offfreq;

            figure; 
            plot(1:31, lvl_drnl(j,:),'LineWidth',2), grid on, hold on
            plot(1:31, lvl_gt(j,:),'r')
            title(sprintf('Excitation pattern in DRNL filterbank (no outer- nor middle - ear compensation)\n%s centred at fc = %.1f [Hz]',tonelabel,fton(j)));
            legend('DRNL','4th-order gt','Location','SouthEast')

            xlabel('Frequency [Hz]')
            ylabel('Log-spectrum at DRNL output [dB]')

            set(gca,'XTick',[1:4:31]);
            set(gca,'XTickLabel',round(audtofreq([1:4:31]+2)));

            ha2(end+1) = gca;
            hFig3(end+1) = gcf;
        end  

    end

    linkaxes(ha2,'xy');
    switch typetone
        case 1
            ylim([SPL-95 SPL+5])
        case 2
            ylim([SPL-95 SPL+15])
        otherwise
            ylim([SPL-95 SPL+5])
    end

    xlim([0.5 31+0.5])

    if bSave
        for i = 1:length(hFig3) 
            Saveas(hFig3(i),sprintf('%sfig-DRNL-no-om-TF-%.0f-SPL-%.0f-dB-%s',dirout,i,SPL,tonelabel_short));
        end
        % close all
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the report:
    % 1. With outer- and middle- ear TFs:
    % fc, Total level, level linear part, level non-linear part, contribution of non-linear path
    % different rows are different centre frequencies
    % var2latex( round(10*[ft' lvldrnl])/10 )

    % 2. Without outer- and middle- ear TFs:
    % fc, Total level, level linear part, level non-linear part, contribution of non-linear path
    % var2latex( round(10*[ft' lvldrnl_noear])/10 );

    % % 3. fc, When introducing TFs:
    % var2latex( round(10*[ft' lvldrnl_noear(:,1:3)-lvldrnl(:,1:3) lvldrnl(:,4)-lvldrnl_noear(:,4)])/10 );

    % Mixing 1, 2 and 3:
    %                                              % increase in total level and in non-linearity compared to no-ear
    var2latex( round(10*[fton' lvltotdrnl lvldrnl_noear(:,[1 4]) lvltotdrnl(:,[1 4])-lvldrnl_noear(:,[1 4])])/10 );

for j = 1:step_bands:nfc
        if typetone == 1
            insig = insigS(:,j);
        elseif typetone == 2
            insig = insigBW20(:,j);
        elseif typetone == 3
            insig = insigBW100(:,j);
        end
        % Each row is the total output of each GN centred at different fcs:
        WV(j,1) = Get_envelope_descriptors(insig,'W');
        WV(j,2) = Get_envelope_descriptors(insig,'V');
        WV_drnl(j,1) = Get_envelope_descriptors( out_drnl_t(:,j),'W');
        WV_drnl(j,2) = Get_envelope_descriptors( out_drnl_t(:,j),'V');
        WV_drnl_noear(j,1) = Get_envelope_descriptors( out_drnl_noear_t(:,j), 'W');
        WV_drnl_noear(j,2) = Get_envelope_descriptors( out_drnl_noear_t(:,j), 'V');
        WV_gt(j,1) = Get_envelope_descriptors( out_gt_t(:,j) ,'W');
        WV_gt(j,2) = Get_envelope_descriptors( out_gt_t(:,j) ,'V');
        Vtot_on_off(1,1) = Get_envelope_descriptors( out_drnl_t    ,'V' );
        Vtot_on_off(2,1) = Get_envelope_descriptors( out_drnl_on_t ,'V' );
        Vtot_on_off(3,1) = Get_envelope_descriptors( out_drnl_off_t,'V' );
        Vtot_on_off_gt(1,1) = Get_envelope_descriptors( out_gt_t    ,'V' );
        Vtot_on_off_gt(2,1) = Get_envelope_descriptors( out_gt_on_t ,'V' );
        Vtot_on_off_gt(3,1) = Get_envelope_descriptors( out_gt_off_t,'V' );

    end  
    Vtot = [WV(:,2) WV_drnl(:,2) WV_drnl_noear(:,2)  WV_gt(:,2)];

    % V closer to 0 is more fluctuating, approaching -infinity fluctuate les
    
    PVdiffs = [WV_drnl(:,2)-WV(:,2)  WV_drnl_noear(:,2)-WV(:,2)  WV_gt(:,2)-WV(:,2)];

    % at 1300-Hz  on frequency
    var2display = [ft2plot(idxon) lvl_drnl(idxon) ; ft2plot(idxoff) lvl_drnl(idxoff); NaN lvltotdrnl(1)];
    var2display = [var2display Vtot_on_off Vtot_on_off_gt];
    var2latex( [round(var2display(:,1:2)*10)/10 round(var2display(:,3:end)*100)/100] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bAnalysisTemplates == 1
    
    [xx mfc] = modfilterbank_debug([1 0],44100,2000); 
    
	[hFig data] = exp_heijden1995;
    
    bListSupra = 1; % this was done to know which where the suprathreshold levels, 10 dB above
    if bListSupra
        disp('Supra levels: 60,72,84 dB')
        supradiff = 10;
        levSupra = data.Thresholds + supradiff;
        
        for i = 1:5
            fprintf('%s: (%.0f, %.0f, %.0f)\n',data.signaltypes{i},levSupra(i,1),levSupra(i,3),levSupra(i,5));
        end
    end
    
    % fnames = {  '84','m04','m04','m04'; ... % for values computed on 02/09/2015
    %             '72','m22','m22','m22'; ...
    %             '60','m30','m30','m30'};
    
    %                 T , 100-GN,20-GN,100-MN,20-MN
    fnames = {  '84','m08','m21','m25','m29','m34'; ...
                '72','m26','m35','m38','m39','m41'; ...
                '60','m36','m38','m38','m38','m38'};
    iTimes = size(fnames,1);
    h1 = [];
    
    for i = 1:iTimes
    
        idxtone = 3;
        exp1 = sprintf('tmp = load(''template-GN-100-Hz-%s-dB-supra-%s-at-2000-Hz.mat'');',fnames{i,1},fnames{i,idxtone}); % tmp = load('template-GN-100-Hz-84-dB-supra-m04-at-2000-Hz.mat');
        eval(exp1);
        templateGN100 = tmp.template;
        
        idxtone = 4;
        exp1 = sprintf('tmp = load(''template-GN-020-Hz-%s-dB-supra-%s-at-2000-Hz.mat'');',fnames{i,1},fnames{i,idxtone}); % tmp = load('template-GN-020-Hz-84-dB-supra-m04-at-2000-Hz.mat');
        eval(exp1);
        templateGN020 = tmp.template;
        
        idxtone = 2;
        exp1 = sprintf('tmp = load(''template-S-%s-dB-supra-%s-at-2000-Hz.mat'');',fnames{i,1},fnames{i,idxtone}); % tmp = load('template-S-84-dB-supra-m04-at-2000-Hz.mat');
        eval(exp1);
        templateS  = tmp.template;
        
        idxtone = 5;
        exp1 = sprintf('tmp = load(''template-MN-100-Hz-%s-dB-supra-%s-at-2000-Hz.mat'');',fnames{i,1},fnames{i,idxtone}); % tmp = load('template-GN-100-Hz-84-dB-supra-m04-at-2000-Hz.mat');
        eval(exp1);
        templateMN100 = tmp.template;
        
        idxtone = 6;
        exp1 = sprintf('tmp = load(''template-MN-020-Hz-%s-dB-supra-%s-at-2000-Hz.mat'');',fnames{i,1},fnames{i,idxtone}); % tmp = load('template-GN-020-Hz-84-dB-supra-m04-at-2000-Hz.mat');
        eval(exp1);
        templateMN020 = tmp.template;
        
        M = length(mfc);
        N = length(templateGN100)/M;

        t = (1:length(templateGN100))/fs;
        % figure;
        % plot(   t, templateS, ...
        %         t, templateGN100, ...
        %         t, templateGN020, ...
        %         t, templateMN100, ...
        %         t, templateMN020 )
        % legend('T','100-Hz GN','20-Hz GN','100-Hz MN','20-Hz MN')
        figure;
        subplot(1,2,1)
        plot(   t, templateS, ...
                t, templateGN100, ...
                t, templateGN020 )
        legend('T','100-Hz GN','20-Hz GN')
        
        subplot(1,2,2)
        plot(   t, templateS, ...
                t, templateMN100, ...
                t, templateMN020 )
        legend('T','100-Hz MN','20-Hz MN')
        h1(end+1) = gcf;
        
        grid on
        templateGN100 = reshape(templateGN100,N,M);
        templateGN020 = reshape(templateGN020,N,M);
        templateS = reshape(templateS,N,M);
        templateMN100 = reshape(templateMN100,N,M);
        templateMN020 = reshape(templateMN020,N,M);
        
        pGN100(i,:) = sum(templateGN100.*templateGN100)*100;
        pGN020(i,:) = sum(templateGN020.*templateGN020)*100;
        pS(i,:)     = sum(templateS.*templateS)*100;
        pMN100(i,:) = sum(templateMN100.*templateMN100)*100;
        pMN020(i,:) = sum(templateMN020.*templateMN020)*100;

    end
    
    plotOpts.bAddVertical = 0;
    hh = Figure2paperfigureT2(h1,3,2,plotOpts);
    
    
    idxlvlSupra = [5 3 1];
    ha = [];
    figure;
    for i = 1:iTimes
        subplot(3,1,i)
        plot(pS(i,:),'bs-'), hold on, grid on
        plot(pGN100(i,:),'mo-.','LineWidth',1.5)
        plot(pGN020(i,:),'ro-','LineWidth',2,'MarkerFaceColor','r')
        plot(pMN100(i,:),'k>-')
        plot(pMN020(i,:),'k<-','LineWidth',2,'MarkerFaceColor','k')
        
        title(sprintf('Masker level = %s dB SPL, test level at suprathreshold: (%.0f,%.0f,%.0f,%.0f,%.0f) dB',fnames{i,1},levSupra(:,idxlvlSupra(i)))) %,levSupra(i,3),levSupra(i,5)))
        if i == 1
            legend(data.signaltypes,'Location','NorthEast') % legend('T','100-Hz GN','20-Hz GN')
        end
        ha(end+1) = gca;
        ylabel('Energy [\%]')

    end
    linkaxes(ha,'xy');
    ylim([0 45])
    xlim([0.5 10.5])
    set(ha,'XTick',1:10);
    set(ha,'XTickLabel',round(mfc));

    xlabel('Centre frequency of the band f_c [Hz]')
    
        
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
