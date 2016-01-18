function r20160113_lilliput_example
% function r20160113_lilliput_example
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 13/01/2016
% Last update on: 13/01/2016 
% Last use on   : 13/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FontSize = 14;
close all
dir = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20160125-27-BATWOMAN-Eindhoven\Sounds\';

bPlots4presentation = 0;
bImprovingICRA = 1;
bPlaySounds = 0;

if bPlots4presentation

    snr        = [-12 -8 -4 0 2:7 10];
    snr_idx    = 1:length(snr);
    
    % snrLIST    = [0   2  3  4  5  6  7];
    scoresLIST = [NaN NaN NaN 10 30 30 50 50 60 60 NaN]; % Jan Leys
    scoresHH_t = nan(size(snr));
    scoresAS_t = nan(size(snr));
    idxHH = [3 11];
    idxAS = [1 2 11];
    
    scores_total = [11	14	 8	28	15;
                    15	15	15	30	15];
	scores_total = scores_total(1,:)./scores_total(2,:);
    scoresHH_t(idxHH) = scores_total(1:2)*100;
    scoresAS_t(idxAS) = scores_total(3:5)*100;
    
    snrAS      = [-12  -8  10];
    scoresAS   = [ 30  90 100; ... % i - ie
                  100 100 100];
    snrHH      = [-4  10];
    scoresHH   = [ 60 100; ... % i - ie
                  100 100];
              
	figure;
    plot(snr_idx,scoresLIST,'ro-','LineWidth',2); grid on
    Xlabel('SNR [dB]',FontSize)
    Ylabel('Scores [%]',FontSize)
    hold on;
    ylim([0 105])
    disp('')
    
    plot(snr_idx,scoresAS_t,'ks','LineWidth',2);
    % plot(snr_idx,scoresHH_t,'b>','LineWidth',2);
    
    figure;
    plot(scoresAS_t([1:2 11]),'ks-','LineWidth',2); grid on; hold on
    % plot(snr_idx,scoresHH_t,'b>','LineWidth',2);
    legend('Total score','scores: ''i-ie''','scores: ''oo''','Location','SouthEast')
    xlim([0.5 3.5])
    ylim([25 105])
    Xlabel('SNR [dB]',FontSize)
    Ylabel('Scores [%]',FontSize)
    Title('Results for Participant 1')
    set(gca,'XTick',[1 2 3])
    set(gca,'XTickLabel',[-12 -8 10])
    set(gca,'FontSize',FontSize)
    hFig = gcf;
    ofile = 'CVC-scores';
    Save_figure_as(hFig,ofile,'emf');
    Save_figure_as(hFig,ofile,'fig');
    
    plot(scoresAS(1,:),'b>--','LineWidth',2); 
    plot(scoresAS(2,:),'r<-.','LineWidth',2); 
    legend('Total score','scores: ''i-ie''','scores: ''oo''','Location','SouthEast')
        
    hFig = gcf;
    Save_figure_as(hFig,[ofile '-detailed'],'emf');
    Save_figure_as(hFig,[ofile '-detailed'],'fig');
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bImprovingICRA
    
    % Sounds to be processed:
    note_test = 'A4';
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];
    
    fname1suffix = 'GRAF28-A4_3';
    fname1   = [dir_where fname1suffix '.wav'];
    [insig fs] = Wavread(fname1);
    [icra info4noise] = icra5_noise4piano(insig,fs);
    
    t = ( 1:length(insig))/fs;
    %%%
    % % Uncomment the following lines to generate a 5-sec noise buffer:
    noisebuf   = AM_random_noise(20,8000,60,5,fs);
    pede(:,1) = il_get_othe_noise(noisebuf, info4noise);
    pede(:,2) = il_get_othe_noise(noisebuf, info4noise);
    pede(:,3) = il_get_othe_noise(noisebuf, info4noise);
    
    icra_pede(:,1) = icra+From_dB(-10)*pede(:,1);
    icra_pede_e(:,1) = il_get_envelope(icra_pede(:,1),fs);
    
    ha = [];
    figure;
    subplot(4,1,1)
    plot(t,icra_pede); hold on, grid on
    plot(t, icra_pede_e,'k','LineWidth',2);
    plot(t,-icra_pede_e,'k','LineWidth',2);
    ha(end+1) = gca;
    
    icra_pede(:,2) = icra+From_dB(-20)*pede(:,2);
    icra_pede_e(:,2) = il_get_envelope(icra_pede(:,2),fs);
    subplot(4,1,2)
    plot(t, icra_pede(:,2)); hold on, grid on
    plot(t, icra_pede_e(:,2),'k','LineWidth',2);
    plot(t,-icra_pede_e(:,2),'k','LineWidth',2);
    ha(end+1) = gca;
    
    icra_pede(:,3) = icra+From_dB(-30)*pede(:,3);
    icra_pede_e(:,3) = il_get_envelope(icra_pede(:,3),fs);
    subplot(4,1,3)
    plot(t, icra_pede(:,3)); hold on, grid on
    plot(t, icra_pede_e(:,3),'k','LineWidth',2);
    plot(t,-icra_pede_e(:,3),'k','LineWidth',2);
    ha(end+1) = gca;
    
    icra_pede_e(:,4) = il_get_envelope(icra,fs);
    subplot(4,1,4)
    plot(t, icra); hold on, grid on
    plot(t, icra_pede_e(:,4),'k','LineWidth',2);
    plot(t,-icra_pede_e(:,4),'k','LineWidth',2);
    ha(end+1) = gca;
        
    linkaxes(ha,'y');
    ylim([-0.2 0.2])
    
    figure;
    plot(t,To_dB(icra_pede_e)); hold on 
    
    MaxdB = max(To_dB(icra_pede_e(:,end)));
    LimsT = minmax(t);
    lines2plot = [MaxdB-10  MaxdB-10; ...
                  MaxdB-20  MaxdB-20; ...
                  MaxdB-30  MaxdB-30];
    plot(LimsT,lines2plot); 
    
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bPlaySounds

% filen = 'weightednoise_nosilence28.4.wav';
filen = 'bruin.wav';
file1 = 'wdpiet.wav';
file2 = 'wdpit.wav';
file3 = 'wdpoot.wav';

[n fs]  = Wavread([dir filen]);
x1      = Wavread([dir file1]);
x2      = Wavread([dir file2]);
x3      = Wavread([dir file3]);

n = n(:,1);

N1 = length(x1);
N2 = length(x2);

sildur = 0.4;
L = 3*sildur*fs + N1 + N2;
sigdur = L/fs;

sil = Gen_silence(sildur,fs);
SNRstart = 10;

while SNRstart >= -5
    
    noi = Randomise_insig(n);
    noi = noi(1:L);
    
    outsig = From_dB(SNRstart)*noi+[sil; x1; sil; x2; sil];
    
    sound(outsig,fs);
    disp(['Press anything to continue...' num2str(-1*SNRstart) ' dB'])
    pause()
    
    SNRstart = SNRstart-3;
        
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function yenv = il_get_envelope(insig,fs,fc)
% This inline function was also used in r20151223_update_Antoine.m

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);

function noise = il_get_othe_noise(noisebuf, info4noise)

noisebuf    = Randomise_insig(noisebuf);
noisebands  = auditoryfilterbank(noisebuf,info4noise.fs);
for i = 1:size(noisebands,2)
    noisebands(:,i) = setdbspl(noisebands(:,i),info4noise.RMS(i)+100);
end
noise = sum(noisebands,2);
noise = setdbspl(noise,info4noise.RMSmax+100);
noise = noise(1:info4noise.Length);
