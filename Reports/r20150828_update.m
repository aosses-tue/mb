function r20150828_update(typetone,SPL)
% function r20150828_update(typetone,SPL)
%
% 1. Description:
%       Computing excitation patterns of the DRNL
%       TODO: 
%           - Middleear Jepsen
%           - Analyse all the bands
%           
% 2. Stand-alone example:
%       typetone = 1; % sine tones
%       SPL = 60; % dB SPL
%       r20150828_update(typetone,SPL)
% 
%       typetone = 2; % 20-Hz Gaussian noise
%       SPL = 60; % dB SPL
%       r20150828_update(typetone,SPL)
%
%       typetone = 3; % 100-Hz Gaussian noise
%       SPL = 60; % dB SPL
%       r20150828_update(typetone,SPL)
%
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 27/08/2015
% Last update on: 27/08/2015 
% Last use on   : 27/08/2015 
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

dirout = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-08-21-update\Figures-eps\';

dBFS = 100;
do_outermiddleear = 1;
fs = 44100;
ft = 3:33; % ERB
ft = audtofreq(ft,'erb');

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
    
    set(gca,'XTick',ft(1:3:end));
    set(gca,'XTickLabel',round(ft(1:3:end)));
    xlim([min(ft) max(ft)])
    
    ylabel('Gain [dB]')
    ylim([-35 28])
    xlabel('Frequency [Hz]')
    
    hFig1(end+1) = gcf;
    
    if bSave
        Saveas(hFig1(end),[dirout 'om-ear']);
    end
end

dur = 1; % s

nfc = length(ft);
nSPL = length(SPL);
lvl = nan(nfc,nfc);
lvltot = nan(nfc,1);

lvl_gt = nan(nfc,nfc);
lvltot_gt = nan(nfc,1);
step_bands = 1; % set to 1 to analyse every band

insigBW100 = [];
insigBW20 = [];
ha1 = [];
ha2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:step_bands:nfc
    
    switch typetone
        case 1 % sinusoid
            insig = Create_sin(ft(j),dur,fs);
            tonelabel = 'Test tone';
            tonelabel_short = 'T';
        case 2 % 20-Hz Gaussian
            BW = 20;
            insig = AM_random_noise_BW(ft(j),BW,60,dur,fs);
            insigBW20 = [insigBW20 insig];
            tonelabel = '20-Hz GN signal';
            tonelabel_short = 'GN-020-Hz';
            W(j) = Get_envelope_descriptors(insig,'W');
            V(j) = Get_envelope_descriptors(insig,'V');
        case 3 % 100-Hz Gaussian
            BW = 100;
            insig = AM_random_noise_BW(ft(j),BW,60,dur,fs);
            insigBW100 = [insigBW100 insig];
            tonelabel = '100-Hz GN signal';
            tonelabel_short = 'GN-100-Hz';
            W(j) = Get_envelope_descriptors(insig,'W');
            V(j) = Get_envelope_descriptors(insig,'V');
    end
    
    for k = 1:length(SPL) 

        insigM = setdbspl( insig,SPL(k) );

        [out_drnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs); 
        lvl(j,:) = rmsdb(out_drnl) + dBFS;
        
        out_drnl_t(:,j) = sum(out_drnl,2);
                
        lvldrnl(j,1) = sum_db( rmsdb(out_drnl) ) + dBFS;
        lvldrnl(j,2) = sum_db( rmsdb(paramsouts.outsiglin) ) + dBFS;
        lvldrnl(j,3) = sum_db( rmsdb(paramsouts.outsignlin) ) + dBFS;
        lvldrnl(j,4) = 100*sum(rms(paramsouts.outsignlin)) / sum(rms(out_drnl));
        
        [out_gt fcgamma]= auditoryfilterbank(insigM,fs);
        lvl_gt(j,:)     = rmsdb(out_gt) + dBFS;
        lvltot_gt(j,1)  = sum_db( rmsdb(out_gt) ) + dBFS;
        
        out_gt_t(:,j) = sum(out_gt,2);
        
        figure; 
        plot(1:31, lvl(j,:),'LineWidth',2), grid on, hold on
        plot(1:31, lvl_gt(j,:),'r')
        title(sprintf('Excitation pattern in DRNL filterbank (outer- and middle - ear compensation)\n%s centred at fc = %.1f [Hz]',tonelabel,ft(j)));
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
        Saveas(hFig2(i),sprintf('%sfig-DRNL-with-om-TF-%.0f-SPL-%.0f-dB-%s',dirout,i,SPL,tonelabel_short)); 
    end
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:step_bands:nfc
    
    switch typetone
        case 1 % sinusoid
            insig = Create_sin(ft(j),dur,fs);
        case 2 % 20-Hz Gaussian
            BW = 20;
            insig = insigBW20(:,j);
        case 3 % 100-Hz Gaussian
            BW = 100;
            insig = insigBW100(:,j);
    end
    
    for k = 1:length(SPL) 

        insigM = setdbspl( insig,SPL(k) );

        % [out_drnl fcdrnl paramsouts] = drnl_CASP_debug(gaindb(insigM,calfactor),fs,'nomiddleear','noouterear');
        [out_drnl fcdrnl paramsouts] = drnl_CASP_debug( insigM,fs,'nomiddleear','noouterear');
        lvl_noear(j,:) = rmsdb(out_drnl) + dBFS;

        out_drnl_noear_t(:,j) = sum(out_drnl,2);
        
        lvldrnl_noear(j,1) = sum_db( rmsdb(out_drnl) ) + dBFS;
        lvldrnl_noear(j,2) = sum_db( rmsdb(paramsouts.outsiglin) ) + dBFS;
        lvldrnl_noear(j,3) = sum_db( rmsdb(paramsouts.outsignlin) ) + dBFS;
        lvldrnl_noear(j,4) = 100*sum(rms(paramsouts.outsignlin)) / sum(rms(out_drnl));
        
        figure; 
        plot(1:31, lvl_noear(j,:),'LineWidth',2), grid on, hold on
        plot(1:31, lvl_gt(j,:),'r')
        title(sprintf('Excitation pattern in DRNL filterbank (no outer- nor middle- ear compensation)\n%s centred at fc = %.1f [Hz]',tonelabel,ft(j)));
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
    close all
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
var2latex( round(10*[ft' lvldrnl lvldrnl_noear(:,[1 4]) lvldrnl(:,[1 4])-lvldrnl_noear(:,[1 4])])/10 );

if (typetone == 2) | (typetone == 3)
    for j = 1:step_bands:nfc
        if typetone == 2
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
        
    end  
end

Vtot = [WV(:,2) WV_drnl(:,2) WV_drnl_noear(:,2)  WV_gt(:,2)];
% V closer to 0 is more fluctuating, approaching -infinity fluctuate les
pV = [prctile(WV(:,2), 5)  prctile(WV_drnl(:,2), 5)  prctile(WV_drnl_noear(:,2), 5)  prctile(WV_gt(:,2), 5); ...
     prctile(WV(:,2),50)  prctile(WV_drnl(:,2),50)  prctile(WV_drnl_noear(:,2),50)  prctile(WV_gt(:,2),50); ...
     prctile(WV(:,2),95)  prctile(WV_drnl(:,2),95)  prctile(WV_drnl_noear(:,2),95)  prctile(WV_gt(:,2),95)];

PVdiffs = [WV_drnl(:,2)-WV(:,2)  WV_drnl_noear(:,2)-WV(:,2)  WV_gt(:,2)-WV(:,2)];
P2 = [prctile(PVdiffs,5); ... 
      prctile(PVdiffs,50); ...
      prctile(PVdiffs,95)]; 

var2latex([round(10*[ft nan(1,3)]')/10 round(1000*[Vtot PVdiffs; pV P2])/1000])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
