function r20160502_figures4seminar
% function r20160502_figures4seminar
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 30/04/2016
% Last update on: 30/04/2016 
% Last use on   : 30/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

dirmain   = ['D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20160503-Acoustics-TUe' delim];
dirthis   = ['Figures-MATLAB'      delim];
dirres    = [dirmain 'Experiments' delim];
diraudio  = [dirmain 'Sounds'      delim];

outputdir = [dirmain dirthis];
Mkdir(outputdir);
Mkdir(diraudio);

FontSize = 14;
bDoICRAspeech = 1;
bFig2 = 0; % Slide 3
bFig3 = 0; % Slide 5 / SRTs for three pairs
bFig5 = 1; % Slide 7 / ICRA + Spectra, change manually the piano you want to plot
% bFig6 = 0;
bFig7 = 0;

bCheckAPEX_sounds = 0; % creates an APEX experiment with two piano sounds to check manually the SNRs

h = [];

if bCheckAPEX_sounds
    dirpi = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Piano\Pilot-ICRA-v2\Stage-9-Pilot-LB-XK\APEX-18-sone\Stimuli-pilot-20160502-B\';
    piano1nonote = 'P0t2';
    piano2nonote = 'P4t3';
    note = 'C4';
    piano1 = [note piano1nonote];
    piano2 = [note piano2nonote];

    p.piano1 = piano1;
    p.piano2 = piano2;
    p.piano1nonote = piano1nonote;
    p.piano2nonote = piano2nonote;
    p.note = note;

    template = [dirpi 'test-XML-red-template.xml.apx302'];
    outputfile = sprintf('%stest-XML-%s-%s.apx302',dirpi,piano1,piano2);
    out2print = readfile_replace(template,p);
    fid     = fopen(outputfile, 'w');
    fwrite(fid, out2print);
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDoICRAspeech
    
    file = [diraudio 'words-seminar.wav'];
    fileo= [diraudio 'words-seminar-ICRA.wav'];
    [insig fs] = Wavread(file);
    method = 1;
    effort_type = 'normal';
    outsig = icra5_noise(insig,fs,method,effort_type);
    lvl = rmsdb(insig);
    outsig = setdbspl(outsig,lvl+100);
    sound(insig,fs);
    pause(length(insig)/fs*1.2);
    sound(outsig,fs);

    Wavwrite(outsig,fs,fileo);
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2 (slide 3)
if bFig2
    h = r20150727_get_simulated_sounds;
    name = 'simsounds';
    Saveas(h(2),sprintf('%s%s',outputdir,name),'epsc')
    Saveas(h(2),sprintf('%s%s',outputdir,name),'emf')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bFig3
    reve = [];
    thres = [];
    N4ave = 8;

    files = {'piano_multi-PILOT-Pres-C4-adaptive-AO.apr', ...
             'piano_multi-PILOT-Pres-C4-adaptive-AO-1.apr', ...
             'piano_PILOT-Pres-C4-adaptive-02-AO.apr'};
    
    for i = 1:length(files)
        [SRT xx xx out] = quick_staircases([dirres files{i}]);
        proc = fieldnames(out.staircases);

        for k = 1:length(proc)
            [reve(end+1,:) xx reve2plot staircase2plot] = Get_mAFC_reversals(out.staircases.(proc{k}));

            reve2use = reve(end,end-N4ave+1:end);
            thres(end+1,1) = median(reve2use);
            thres(end  ,2) = str2num(proc{k}(end-1:end)); % last two characters
            thres(end  ,3) = i;
            thres(end  ,4) = min(reve2use);
            thres(end  ,5) = max(reve2use);
          
        end
    
    end
    
    labels = {'0-2','2-4','4-0'};
    thres(1,:) = thres(7,:);
    thres(7,:) = [];
    
    [xx idx] = sort(thres(:,2));
    thres = thres(idx,:);
    
    exp1 = transpose( thres(1:2:end,1) );
    exp2 = transpose( thres(2:2:end,1) );
    expmean = mean([exp1; exp2]);
    maxExp = max([exp1; exp2]);
    minExp = min([exp1; exp2]);
    errorU = maxExp - expmean;
    errorL = expmean - minExp;
    offsetx = 0.1;
    figure;
    errorbar([1:3],expmean,errorL,errorU,'o--','LineWidth',2), hold on, grid on
    Ylabel('Threshold, SNR [dB]',FontSize)
    Xlabel('Piano pair')
    ha = gca;
    set(ha,'XTick',[1:3]);
    set(ha,'XTickLabel',labels);
    ylim([-9 7])
    set(ha,'FontSize',FontSize);
    h(end+1) = gcf;
    
    name = 'results-adapt-AO';
    Saveas(h(end),sprintf('%s%s',outputdir,name),'epsc')
    Saveas(h(end),sprintf('%s%s',outputdir,name),'emf')
    
end

if bFig5
        
    dur = 1.5;
    dir_out   = [Get_TUe_data_paths('lx_Text') 'lx2015-12-21-update-ICRA-Antoine' delim 'piano-sounds-new' delim];
    Mkdir(dir_out);
    
    %%%
    dir = [Get_TUe_paths('ICRA_Tobias')];
    addpath([dir delim 'Tools'])
    %%%
    
    % Sounds to be processed:
    note_test = 'C4';
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];

    fname0suffix = 'GH05-C4_2';
    fname1suffix = 'JBS36-C4_2';
    fname2suffix = 'JBS51-4544-C4_3';
    
    fname0   = [dir_where fname0suffix '.wav'];
    fname1   = [dir_where fname1suffix '.wav'];
    fname2   = [dir_where fname2suffix '.wav'];
 
    legend0txt = 'GH05 (C4)';
    legend0txtn = 'ICRA noise from GH05 (C4)';
    colour0 = 'm';

    legend1txt = 'JBS36 (C4)';
    legend1txtn = 'ICRA noise from JBS36 (C4)';
    colour1 = 'b';
    
    legend2txt = 'JBS51-4544 (C4)';
    legend2txtn = 'ICRA noise from JBS51-4544';
    colour2 = 'r';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Reading the wav files:
    [signal0 fs] = Wavread(fname0); 
    [signal1 fs] = Wavread(fname1); 
    [signal2 fs] = Wavread(fname2); 

    [L idxL] = min([length(signal0) length(signal1) length(signal2)]);
    signal0 = Do_truncate(signal0,L);
    signal1 = Do_truncate(signal1,L);
    signal2 = Do_truncate(signal2,L);
 
    dur_ramp = 150; % ms
    signal0 = Do_cos_ramp(signal0,fs,0,dur_ramp);
    signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp);
    signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp);
 
    L = round(dur*fs);
    signal0 = signal0(1:L);
    signal1 = signal1(1:L);
    signal2 = signal2(1:L);
    
    t = (1:L) / fs;

    % 2.2. Storing the aligned waveforms:
    fname_out = sprintf('%s%s-demo.wav',diraudio,fname0suffix);
    Wavwrite(signal0,fs,fname_out);
    fname_out = sprintf('%s%s-demo.wav',diraudio,fname1suffix);
    Wavwrite(signal1,fs,fname_out);
    fname_out = sprintf('%s%s-demo.wav',diraudio,fname2suffix);
    Wavwrite(signal2,fs,fname_out);
 
    % method    = 3; % i = 3 = ERB;
    fname_out0 = sprintf('%snoise-%s-ICRA.wav',diraudio,fname0suffix);
    fname_out1 = sprintf('%snoise-%s-ICRA.wav',diraudio,fname1suffix);
    fname_out2 = sprintf('%snoise-%s-ICRA.wav',diraudio,fname2suffix);

    try
        noise0  = Wavread(fname_out0);
        noise1  = Wavread(fname_out1);
        noise2  = Wavread(fname_out2);
    catch
        disp('At least one of the noises was not found...')
        noise0   = icra_noise4piano(signal0,fs);
        Wavwrite(noise0,fs,fname_out0);
        noise1   = icra_noise4piano(signal1,fs);
        Wavwrite(noise1,fs,fname_out1);
        noise2   = icra_noise4piano(signal2,fs);
        Wavwrite(noise2,fs,fname_out2);
    end

    nSig2plot = 2;
    switch nSig2plot
        case 0
            signal = signal0;
            noise  = noise0;
            fnamesuffix = fname0suffix;
            ylimits = [-50 9];
            legendtxt = legend0txt;
            legendtxtn = legend0txtn;
            colour = colour0;
            name1 = 'waveform-GH05';
            name2 = 'spectrum-GH05';
    
        case 1
            signal = signal1;
            noise  = noise1;
            fnamesuffix = fname1suffix;
            ylimits = [-50 9];
            legendtxt = legend1txt;
            legendtxtn = legend1txtn;
            colour = colour1;
            name1 = 'waveform-JBS36';
            name2 = 'spectrum-JBS36';
    
        case 2
            signal = signal2;
            noise  = noise2;
            fnamesuffix = fname2suffix;
            ylimits = [-50 9];
            legendtxt = legend2txt;
            legendtxtn = legend2txtn;
            colour = colour2;
            name1 = 'waveform-JBS51-4544';
            name2 = 'spectrum-JBS51-4544';
    end
    fc = 20;
    yenv  = il_get_envelope(signal,fs,fc);
    yenvn = il_get_envelope(noise,fs,fc);
 
    figure;
    Text_ylabel = 'Amplitude';
    Text_xlabel = 'Time [s]';
    
    subplot(2,1,1)
    plot(t,signal,colour); hold on, grid on
    plot(t,  yenv,'k','LineWidth',2);
    plot(t, -yenv,'k','LineWidth',2);
    ylabel(Text_ylabel)
    title(legendtxt)
    ha = gca;

    subplot(2,1,2)
    plot(t, noise,colour); hold on, grid on
    plot(t,  yenvn,'k','LineWidth',2);
    plot(t, -yenvn,'k','LineWidth',2);
    ylabel(Text_ylabel)
    title(legendtxtn)
    ha(end+1) = gca;

    linkaxes(ha,'xy');
    xlabel(Text_xlabel)
    ylim([-0.6 0.6])
    xlim([0 max(t)])
 
    plotOpts = [];
    plotOpts.I_Width = 10;
    plotOpts.FontSize = FontSize;
    
    h(end+1) = gcf;
    h(end) = Figure2paperfigureT(h(end),2,1,plotOpts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [RMS_se t_se] = rmsdb_sec(signal,fs,10e-3);
    [max_0 idx_0] = max(RMS_se);

    Ni = find(t<t_se(idx_0),1,'last');
    Nf = Ni + round(500e-3*fs);
    
    % % Calculate LTAS
    [ltassP, fHz] = calcLTAS(signal,fs);
    [ltassN     ] = calcLTAS(noise,fs);
    
    dBMax = max( max([ltassP ltassN]) );
    ltassP = ltassP-dBMax;
    ltassN = ltassN-dBMax;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    % subplot(1,2,1);
    semilogx(   fHz,ltassP,colour, ... 
                fHz,ltassN,'-k'); hold on, grid on;
    xlabel('Frequency [Hz]')
    ylabel('Relative amplidude [dB]')
    axis tight;
    XTick = [50 100 250 500 1000 2000 4000 8000];
    set(gca,'XTick'     ,XTick);
    set(gca,'XTickLabel',XTick);
    ha0 = gca;    
    title(name2figname(fnamesuffix));

    linkaxes(ha0,'xy');
    ylim(ylimits)
    xlim([200 2048]);

    h(end+1) = gcf;
    h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);

    legend({'signal','noise'},'location','southwest')
   
    Saveas(h(end-1),sprintf('%s%s',outputdir,name1),'epsc')
    Saveas(h(end-1),sprintf('%s%s',outputdir,name1),'emf')
    
    Saveas(h(end),sprintf('%s%s',outputdir,name2),'epsc')
    Saveas(h(end),sprintf('%s%s',outputdir,name2),'emf')
    
    rmpath([dir delim 'Tools'])  
end

plotOpts = [];
plotOpts.I_Width = 12;
plotOpts.FontSize = FontSize;

if bFig7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dir = 'D:\Databases\dir04-Psychoacoustics\WAE\saves-my-results\';
    files = {'save-AO-silence-order.xml','save-AO-silence-1.xml'};
    [score_final score_final_pair] = r20160429_experiments_WAE(dir, files);
    
    label_piano={'0','2t1','2t2','4','6'};
    for i = 1:size(score_final_pair,1)
        p1 = score_final_pair(i,2);
        p2 = score_final_pair(i,3);
        label_pair{i} = sprintf('%s/%s',label_piano{p1},label_piano{p2});
    end
    
    [label_pair idx] = sort(label_pair);
    score_pair = score_final_pair(idx,:);
    [xx idx2] = sort(score_pair(:,1),'descend');
    score_pair = score_pair(idx2,:);
    label_pair = label_pair(idx2);
    
    offsetx = 0.1;
    figure;
    xdata = 1:size(score_pair);
    plot(xdata-offsetx,score_pair(:,1),'ro','LineWidth',2); grid on
    Xlabel('Piano pair',FontSize)
    Ylabel('Most similar pair, Score[\%]',FontSize)
    
    xlim([min(xdata)-0.5 max(xdata)+0.5]);
    ylim([0 103])
    set(gca,'XTick',xdata);
    h(end+1) = gcf;
    h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
    set(gca,'XTickLabel',label_pair);
    set(gca,'FontSize',FontSize);
    hold on;
    legend('silence')
    
    name = 'triadic-one-participant';
    Saveas(h(end),sprintf('%s%s',outputdir,name),'epsc')
    Saveas(h(end),sprintf('%s%s',outputdir,name),'emf')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = {'save-AO-SNR-0-dB.xml','save-AO-SNR-0-dB-1.xml'};
    [score_final score_final_pair] = r20160429_experiments_WAE(dir, files);
    
    [xx idx] = sort(score_final_pair(:,1),'descend');
    
    label_piano={'0','2t1','2t2','4','6'};
    for i = 1:size(score_final_pair,1)
        p1 = score_final_pair(i,2);
        p2 = score_final_pair(i,3);
        label_pair{i} = sprintf('%s/%s',label_piano{p1},label_piano{p2});
    end
    
    [label_pair idx] = sort(label_pair);
    score_pair = score_final_pair(idx,:);
    score_pair = score_pair(idx2,:);
    
    plot(xdata+offsetx,score_pair(:,1),'bs','LineWidth',2); grid on
    
    h(end+1) = gcf;
    % h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
    legend('silence','in noise, SNR = 0 dB')
    name = 'triadic-one-participant-SNR';
    Saveas(h(end),sprintf('%s%s',outputdir,name),'epsc')
    Saveas(h(end),sprintf('%s%s',outputdir,name),'emf')
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);

