function r20160601_piano_figs4paper
% function r20160601_piano_figs4paper
%
% 1. Description:
%       Processing as in r20160502_figures4seminar (for seminar Acoutics at TU/e)
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20160502_figures4seminar
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 21/06/2016
% Last update on: 21/06/2016 
% Last use on   : 21/06/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

dirmain   = ['D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-11-ISMRA' delim];
dirthis   = ['Figures-MATLAB'      delim];
diraudio  = [dirmain 'Sounds'      delim];

outputdir = [dirmain dirthis];
Mkdir(outputdir);
Mkdir(diraudio);

FontSize = 14;
bSave = 0;

bPart1 = 1; % Waveforms          + envelopes (20-Hz LPF)
bPart1b = 1; % Waveforms [dB SPL] + envelopes 
bPart2 = 1; % FFT (LTAS)

h = [];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 2 (slide 3)
% if bFig2
%     h = r20150727_get_simulated_sounds;
%     name = 'simsounds';
%     Saveas(h(2),sprintf('%s%s',outputdir,name),'epsc')
%     Saveas(h(2),sprintf('%s%s',outputdir,name),'emf')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dur = 1.5;
dir_out   = [Get_TUe_data_paths('lx_Text') 'lx2015-12-21-update-ICRA-Antoine' delim 'piano-sounds-new' delim];
Mkdir(dir_out);

% Sounds to be processed:
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '06-Exp-TUe-1-similarity' delim '02-to-be-used-in-AB-comparison' delim];

fname0suffix = 'GH05-Cd5_1-18-sone-dur-1300-ms';
fname1suffix = 'JBS36-Cd5_2-18-sone-dur-1300-ms';
fname2suffix = 'JBS51-4544-Cd5_6-18-sone-dur-1300-ms';

fname0   = [dir_where fname0suffix '.wav'];
fname1   = [dir_where fname1suffix '.wav'];
fname2   = [dir_where fname2suffix '.wav'];

legend0txt = 'GH05 (C#_5)';
legend0txtn = 'ICRA noise from GH05 (C#_5)';
colour0 = 'm';

legend1txt = 'JBS36 (C#_5)';
legend1txtn = 'ICRA noise from JBS36 (C#_5)';
colour1 = 'b';

legend2txt = 'JBS51-4544 (C#_5)';
legend2txtn = 'ICRA noise from JBS51-4544';
colour2 = 'r';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Reading the wav files:
[signal0 fs] = Wavread(fname0); 
[signal1 fs] = Wavread(fname1); 
[signal2 fs] = Wavread(fname2); 

[L idxL] = min([length(signal0) length(signal1) length(signal2)]);


t = (1:L) / fs;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hname = [];

for nSig2plot = 0:2;

    switch nSig2plot
        case 0
            signal = signal0;
            noise  = noise0;
            fnamesuffix = fname0suffix;
            ylimits = [-55 5];
            legendtxt = legend0txt;
            legendtxtn = legend0txtn;
            colour = colour0;
            pianoID = '1';

        case 1
            signal = signal1;
            noise  = noise1;
            fnamesuffix = fname1suffix;
            ylimits = [-55 5];
            legendtxt = legend1txt;
            legendtxtn = legend1txtn;
            colour = colour1;
            pianoID = '3';

        case 2
            signal = signal2;
            noise  = noise2;
            fnamesuffix = fname2suffix;
            ylimits = [-55 5];
            legendtxt = legend2txt;
            legendtxtn = legend2txtn;
            colour = colour2;
            pianoID = '5';
    end
    
    if bPart1
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
        plotOpts.I_Width = 11;
        plotOpts.FontSize = FontSize;

        h(end+1) = gcf;
        hname{end+1} = ['waveform' pianoID];
        h(end) = Figure2paperfigureT(h(end),2,1,plotOpts);
        if bSave
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'epsc');
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'emf');
        end
    end
    
    if bPart1b
        
        % N = 1000;
        % inene = [zeros(N,1); signal; zeros(N,1)];
        inene = signal;
        ene_v1 = il_get_envelope(       signal,fs, 20);
        ene_v2 = il_get_envelope( inene.*inene,fs, 20);
        ene_v3 = il_get_envelope(max(signal,0),fs, 20); % Half-wave rectification
        
        % ene_v2 = ene_v2(N+1:end-N);
        
        plotOpts = [];
        plotOpts.I_Width = 10;
        plotOpts.FontSize = FontSize;
        plotOpts.I_KeepTicks = 0;
        
        figure;
        subplot(1,2,1)
        plot(t,20*log10(abs(2*signal)/2e-5),colour); hold on
        ylabel('Sound Pressure Level [dB]');
        xlabel('Time [s]')
        title(sprintf('piano %s',pianoID));
        
        ylim([30 89]);
        xlim([0 1.3]);
        
        if nSig2plot == 0
            plot(t,20*log10((2*ene_v1)/2e-5),'k','LineWidth',2);
            plot(t,20*log10((2*ene_v2)/2e-5),'k');
            plot(t,20*log10((2*ene_v3)/2e-5),'k--','Linewidth',2);
            legend('Waveform','1. Hilb+LPF','2. Energy+LPF','3. HWR+LPF');
        else
            plot(t,20*log10((2*ene_v1)/2e-5),'k','LineWidth',2);
            legend('Waveform','1. Hilb+LPF');
        end
        
        inene = noise;
        ene_v1 = il_get_envelope(       noise,fs, 20);
        ene_v2 = il_get_envelope( inene.*inene,fs, 20);
        ene_v3 = il_get_envelope(max(noise,0),fs, 20);
        
        %%%
        subplot(1,2,2)
        plot(t,20*log10(abs(2*noise)/2e-5),colour); hold on
        plot(t,20*log10((2*ene_v1)/2e-5),'k','LineWidth',2);
        % legend('Waveform','1. Hilb+LPF');
        title(sprintf('noise for piano %s',pianoID));
        
        ylim([30 89]);
        xlim([0 1.3]);
        xlabel('Time [s]')
        
        h(end+1) = gcf;
        hname{end+1} = ['waveform-SPL-' pianoID];
        h(end) = Figure2paperfigureT(h(end),1,2,plotOpts);
        
        if bSave
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'epsc');
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'emf');
        end
    end
    
    if bPart2

        dir = [Get_TUe_paths('ICRA_Tobias')];
        addpath([dir delim 'Tools'])


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
        % title(name2figname(fnamesuffix));
        title( sprintf('Spectrum: piano %s \nand its ICRA noise',pianoID)  )

        linkaxes(ha0,'xy');
        ylim(ylimits)
        xlim([400 4200]);

        h(end+1) = gcf;
        plotOpts = [];
        plotOpts.I_Width = 11;
        plotOpts.FontSize = FontSize;
        h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
        hname{end+1} = ['spectrum-' pianoID];
        legend({'piano','noise'},'location','northeast')

        if bSave
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'epsc')
            Saveas(h(end),sprintf('%s%s',outputdir,hname{end}),'emf')
        end

        rmpath([dir delim 'Tools'])  

    end
    
    disp('')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);

