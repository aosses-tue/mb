function r20160404_update_proc_ICRA_APEX
% function r20160404_update_proc_ICRA_APEX
%
% 1. Description:
%       Analysing the intermediate results of the ICRA pilot.
% 
% 2. Stand-alone example:
%       r20160404_update_proc_ICRA_APEX;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 01/04/2016
% Last update on: 01/04/2016 
% Last use on   : 01/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDoPlots = 1;
bDoAnaPedestal = 0;
bDoWaveforms = 1;

bSave = 0;

dir_main = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-5-Multiprocedure-auto' delim];
dir_where = [dir_main 'S00-Initial-test-AO' delim];

dirS = {[dir_main 'S00-Initial-test-AO' delim]; ...
        [dir_main 'S01-Initial-test-AS' delim]; ...
        [dir_main 'S02-Initial-test-HW' delim]};

if bDoPlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = {   'piano_multi_result-C2-test-1.apr', 'piano_multi_result-A4-test-1.apr', 'piano_multi_Csh5-test-1.apr'; ...
        'piano_multi_C2-AS.apr', 'piano_multi_A4-AS.apr', 'piano_multi_Csh5-AS.apr'; ...
        'piano_multi_C2-HW.apr', 'piano_multi_A4-HW.apr', 'piano_multi_Csh5-HW.apr'};

h = [];
for j = 1:3
    for i = 1:length(f)
        [SRTtmp stdtmp rest] = quick_staircases([dirS{j} f{j,i}]);
        
        Data17(j,i) = SRTtmp(1); 
        Data71(j,i) = SRTtmp(3);
        Data47(j,i) = SRTtmp(2); 
        Data74(j,i) = SRTtmp(4); 
        
        if bSave
            h(end+1) = Figure2paperfigureT(gcf,2,2);
            Saveas(h(end),sprintf('S%.0f-staircase-%.0f',j-1,i));
        end
    end
end

xoffset = 0.05;
figure;
errorbar((1:3)-2*xoffset,mean(Data17),std(Data17),'ro','LineWidth',2), hold on
errorbar((1:3)+2*xoffset,mean(Data47),std(Data47),'bo' )

errorbar((1:3)-1*xoffset,mean(Data71),std(Data71),'rs','LineWidth',2)
errorbar((1:3)+1*xoffset,mean(Data74),std(Data74),'bs' )

grid on
ha = gca;
set(ha,'XTick',[1 2 3]);
set(ha,'XTickLabel',{'C2','A4','Csh5'});
set(ha,'FontSize',14)
legend('GRAF28','JBS51-4544')
Xlabel('Note',14)
Ylabel('50% SNR (dB)',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bDoAnaPedestal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = [dir_where 'Stimuli-A4' delim 'noise-P1t2_P7t2.wav'];
    fp = [dir_where 'Stimuli-A4' delim 'pedestal-P1t2_P7t2-SNR-10-dB.wav'];

    [xn fs] = Wavread(fn);
    [xp   ] = Wavread(fp);

    xp = From_dB(10)*xp(1:length(xn));

    [dbn t1] = rmsdb_sec(xn,fs,10e-3);
    [dbp t2] = rmsdb_sec(xp,fs,10e-3);

    figure;
    plot(t1,dbn,t2,dbp);
    legend('ICRA','pede')

    [iin,ssn] = Get_Loudness_MGB( xn         , fs);
    [iii,sss] = Get_Loudness_MGB( xp(1:44100), fs);

    figure;
    plot(1:length(ssn),ssn,1:length(sss),sss);
    legend('ICRA','pede')

    disp('')

end

if bDoWaveforms
    
    dire = ['D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-04-01-update_ICRA_APEX\Audio' delim];
    tlims = [0.03 1.13];
    clICRA = rgb('DarkKhaki');
    xlim(tlims)
    
    [x1 fs] = Wavread([dire 'P7t2.wav']);
    t1  = ( 1:length(x1) )/fs;
    
    pr  = 2*x1;
    SPL = 20*log10(abs(pr)/2e-5);
    
    figure;
    subplot(3,1,1)
    plot(t1,pr,'k');
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    grid on;
    ha = gca;
    title('NS19, t2')
    xlim(tlims)
    
    subplot(3,1,2)
    plot(t1,SPL,'k');
    xlabel('Time [s]')
    ylabel('SPL [dB]')
    grid on, hold on;
    xlim(tlims)
    
    y1 = abs(pr); % pr.*pr;
    y2 = Apply_IIR_Butter(y1,fs,20,'low',4); 
    plot(t1,20*log10(2*y2/2e-5),'r','LineWidth',2) % twice because I am deleting negative part
    ylim([33 85])   
    ha(end+1) = gca;
    
    subplot(3,1,3)
    nfft = 8192;
    wlen = nfft/2;
    overlap = 87.5; % percent
    nwtype = 4; % 4 = hamming
    opts.bColourBar = 0;
    stft(x1,fs, nfft, wlen, overlap, nwtype,opts);
    fmin = 400;
    fmax = fmin*11;
    ylim([fmin fmax])
    ha(end+1) = gca;
    title('')
    xlabel('Time [s]')
    xlim(tlims)
    
    linkaxes(ha,'x');
    
    opts = [];
    h = Figure2paperfigureT(gcf,3,1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [xn fs] = Wavread([dire 'noise-P4t4_P7t2.wav']);
    
    prn  = 2*xn;
    SPLn = 20*log10(abs(prn)/2e-5);
    
    figure;
    subplot(3,1,1)
    plot(t1,prn,'Color',clICRA);
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    grid on;
    ha(end+1) = gca;
    title(sprintf('ICRA (pair NS19 t2, JBS51-4544 t4)\n SNR = 0 dB'))
    xlim(tlims)
    
    subplot(3,1,2)
    plot(t1,SPLn,'Color',clICRA);
    xlabel('Time [s]')
    ylabel('SPL [dB]')
    grid on, hold on;
    xlim(tlims)
    
    y1 = abs(prn); % pr.*pr;
    y2 = Apply_IIR_Butter(y1,fs,20,'low',4); 
    plot(t1,20*log10(2*y2/2e-5),'Color',rgb('Gray'),'LineWidth',2) % twice because I am deleting negative part
    ylim([33 85])   
        
    ha(end+1) = gca;
    
    subplot(3,1,3)
    nfft = 8192;
    wlen = nfft/2;
    overlap = 87.5; % percent
    nwtype = 4; % 4 = hamming
    opts.bColourBar = 0;
    stft(xn,fs, nfft, wlen, overlap, nwtype,opts);
    fmin = 400;
    fmax = fmin*11;
    ylim([fmin fmax])
    ha(end+1) = gca;
    title('')
    xlabel('Time [s]')
    xlim(tlims)
    
    % linkaxes(ha,'x');
    
    opts = [];
    h(end+1) = Figure2paperfigureT(gcf,3,1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [xp fs] = Wavread([dire 'pedestal-P4t4_P7t2-SNR-10-dB.wav']);
    
    SNR = 8; 
    SNRpede = 20;
    prp    = From_dB(-SNRpede)*2*xp(1:length(x1));
    prnAtt = From_dB(-SNR)*prn;
    
    
    SPLp = 20*log10(abs(prp)/2e-5);
    SPLnAtt = 20*log10(abs(prnAtt)/2e-5);
    
    figure;
    subplot(3,1,1)
    plot(t1,pr,'k'); hold on
    plot(t1,prnAtt,'Color',clICRA)
    plot(t1,prp,'Color',rgb('YellowGreen'));
    xlabel('Time [s]')
    ylabel('Pressure [Pa]')
    grid on;
    ha(end+1) = gca;
    title(sprintf('Interval NS19 t2 + ICRA + Pedestal \n (SNR = %.0f dB)',SNR))
    xlim(tlims)
    
    subplot(3,1,2)
    plot(t1,SPL,'k'); hold on
    plot(t1,SPLnAtt,'Color',clICRA);
    plot(t1,SPLp,'Color',rgb('YellowGreen'));
    
    xlabel('Time [s]')
    ylabel('SPL [dB]')
    grid on, hold on;
    xlim(tlims)
    
    y1 = abs(prnAtt); % pr.*pr;
    y2 = Apply_IIR_Butter(y1,fs,20,'low',4); 
    plot(t1,20*log10(2*y2/2e-5),'Color',rgb('Gray'),'LineWidth',2) % twice because I am deleting negative part
    
    y1 = abs(pr); % pr.*pr;
    y2 = Apply_IIR_Butter(y1,fs,20,'low',4); 
    plot(t1,20*log10(2*y2/2e-5),'r','LineWidth',2) % twice because I am deleting negative part
    
    y1 = abs(prp); % pr.*pr;
    y2 = Apply_IIR_Butter(y1,fs,20,'low',4); 
    plot(t1,20*log10(2*y2/2e-5),'Color',rgb('LightGreen'),'LineWidth',2) % twice because I am deleting negative part
        
    ylim([33 85])   
        
    ha(end+1) = gca;
    
    subplot(3,1,3)
    nfft = 8192;
    wlen = nfft/2;
    overlap = 87.5; % percent
    nwtype = 4; % 4 = hamming
    opts.bColourBar = 0;
    stft(x1 + From_dB(-SNR)*xn + From_dB(-SNRpede)*xp(1:length(x1)),fs, nfft, wlen, overlap, nwtype,opts);
    fmin = 400;
    fmax = fmin*11;
    ylim([fmin fmax])
    ha(end+1) = gca;
    title('')
    xlabel('Time [s]')
    xlim(tlims)
    
    % linkaxes(ha,'x');
    
    opts = [];
    h(end+1) = Figure2paperfigureT(gcf,3,1,1);
    
    for i = 1:3
        Saveas(h(i),['waveforms-' num2str(i)],'epsc');
    end
    
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


