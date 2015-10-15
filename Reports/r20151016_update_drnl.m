function y = r20151016_update_drnl(bDoParts)
% function y = r20151016_update_drnl(bDoParts)
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

if nargin == 0
    bDoParts = [1 0 0 0 0];
end

bDoTFDRNL2  = bDoParts(1);

h = [];
dire = [Get_TUe_paths('outputs') 'AMTControl-examples' delim]; % 'D:\Output\AMTControl-examples\'
dBFS = 100;

fmax2plot = 1000;
N = 8192*8;
K = N/2;
fs = 44100;

bFig2ab = 1;
bFig2c = 0;
    
if bDoTFDRNL2
    
    SPL = [0 20:10:90 100];
    SPL_fig2c = [30 60 90]; % dB
    bandidx = [1 5 9 14 25 31]; % 250, 500, 1000, 4000 Hz respectively
    
    % fc_fig2a = [100 250  500 1000 4000 8000];
    fc_fig2a = [86.9 257.4 520 1055.8  3982.6 7819.2];
    fc_fig2b = [1000 2400 4000 8000];
    % fc_fig2c = Get_OB_freqs(3,250,2000);
    dur = N/fs;
    
    % %% Generation of Figure 2.c-e
    % if bFig2c
    %     nfc = length(fc_fig2c);
    %     for j = 1:nfc
    %         insig = Create_sin(fc_fig2c(j),dur,fs);
    %         for k = 1:length(SPL_fig2c) 
    % 
    %             insigM = setdbspl( insig,SPL_fig2c(k) );
    % 
    %             [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs);
    %             out = outsigdrnl(:,bandidx(3)); % band centred at 1 kHz
    %             lvl(j,k) = rmsdb(out);
    % 
    %             [outsiggamma fcgamma] = auditoryfilterbank(insigM,fs);
    %             outgamma = outsiggamma(:,bandidx(3)); % band centred at 1 kHz
    %             lvlgamma(j,k) = rmsdb(outgamma);
    %         end       
    % 
    %     end
    %     % Figure 2.c-e
    %     figure; 
    %     subplot(1,3,1)
    %     plot(lvl(:,1)-max(lvl(:,1)),'bo-','LineWidth',2); hold on; 
    %     plot(lvlgamma(:,1)-max(lvlgamma(:,1)),'rx--'); grid on
    %     set(gca,'XTick',[1:2:nfc]);
    %     set(gca,'XTickLabel',round(fc_fig2c(1:2:end)));
    %     xlabel('Frequency [Hz]')
    %     ylabel('DRNL output [dB re max]')
    % 
    %     title('Iso-intensity response functions')
    %     subplot(1,3,2)
    %     plot(lvl(:,2)-max(lvl(:,2)),'bo-','LineWidth',2); hold on; 
    %     plot(lvlgamma(:,2)-max(lvlgamma(:,2)),'rx--'); grid on
    %     set(gca,'XTick',[1:2:nfc]);
    %     set(gca,'XTickLabel',round(fc_fig2c(1:2:end)));
    %     xlabel('Frequency [Hz]')
    %     ylabel('DRNL output [dB re max]')
    % 
    %     subplot(1,3,3)
    %     plot(lvl(:,3)-max(lvl(:,3)),'bo-','LineWidth',2); hold on; 
    %     plot(lvlgamma(:,3)-max(lvlgamma(:,3)),'rx--'); grid on
    %     set(gca,'XTick',[1:2:nfc]);
    %     set(gca,'XTickLabel',round(fc_fig2c(1:2:end)));
    %     xlabel('Frequency [Hz]')
    %     ylabel('DRNL output [dB re max]')
    % 
    %     legend('DRNL','4th-order Gammatone')
    %     h(end+1) = gcf;
    % end
    
    %% Generation of Figure 2.a
    if bFig2ab
        for j = 1:length(fc_fig2a)
            insig = Create_sin(fc_fig2a(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs);    % DRNL
                % [outsiggamma fcgamma]          = auditoryfilterbank(insigM,fs); % Gamma-tone
                % outgamma = outsiggamma(:,bandidx(j)); % band centred at 1 kHz
                
                out = outsigdrnl(:,bandidx(j));
                % lvl2agamma(j,k) = rmsdb(outgamma);
                lvl2a(j,k) = rmsdb(out);
                
                disp('')
            end       

        end

        % Figure 2.a
        figure;
        subplot(1,2,1)
        plot(   SPL,lvl2a(1,:),'k.--','LineWidth',1 ); hold on
        plot(   SPL,lvl2a(2,:),'go-.','LineWidth',2 )
        plot(   SPL,lvl2a(3,:),'r--') % 1k
        plot(   SPL,lvl2a(4,:),'b-', 'LineWidth',2 ); % 4k
        plot(   SPL,lvl2a(5,:),'m>-', 'LineWidth',2 ); % 4k
        legend('250 Hz','500 Hz','1 kHz','4 kHz','8 kHz')
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
        for j = 1:4
            insig = Create_sin(fc_fig2b(j),dur,fs);
            for k = 1:length(SPL) 

                insigM = setdbspl( insig,SPL(k) );

                [outsigdrnl fcdrnl paramsouts] = drnl_CASP_debug(insigM,fs);

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
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Figure 2.a gamma
        % figure;
        % subplot(1,2,1)
        % plot(   SPL,lvl2a(1,:),'k.--','LineWidth',1 ); hold on
        % plot(   SPL,lvl2a(2,:),'go-.','LineWidth',2 )
        % plot(   SPL,lvl2a(3,:),'r--') % 1k
        % plot(   SPL,lvl2a(4,:),'b-', 'LineWidth',2 ); % 4k
        % legend('250 Hz','500 Hz','1 kHz','4 kHz')
        % title('I/O functions, different CFs')
        % xlabel('Input level [dB SPL]')
        % ylabel('DRNL output [dB re 1 m/s]')
        % xlim([0 100])
        % ylim([-100 0])
        % grid on
        % 
        % subplot(1,2,2)
        % plot(   SPL,lvl2agamma(1,:),'k.--','LineWidth',1 ); hold on
        % plot(   SPL,lvl2agamma(2,:),'go-.','LineWidth',2 )
        % plot(   SPL,lvl2agamma(3,:),'r--') % 1k
        % plot(   SPL,lvl2agamma(4,:),'b-', 'LineWidth',2 ); % 4k
        % legend('250 Hz','500 Hz','1 kHz','4 kHz')
        % title('I/O functions, different CFs')
        % xlabel('Input level [dB SPL]')
        % ylabel('Gammatone output [dB]')
        % xlim([0 100])
        % ylim([-100 0])
        % grid on
        h(end+1) = gcf;
    end
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
