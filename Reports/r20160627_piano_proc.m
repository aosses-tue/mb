function r20160627_piano_proc(fname1,fname2)
% function r20160627_piano_proc(fname1,fname2)
%
% 1. Description:
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
bPart2 = 1; % Waveforms [dB SPL] + envelopes 
bPart3 = 1; % FFT (LTAS)

h = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dur = 1.5;
dir_out   = [Get_TUe_data_paths('lx_Text') 'lx2015-12-21-update-ICRA-Antoine' delim 'piano-sounds-new' delim];
Mkdir(dir_out);

if nargin == 0 
    % Sounds to be processed:
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '06-Exp-TUe-1-similarity' delim '02-to-be-used-in-AB-comparison' delim];

    fname1suffix = 'JBS36-Cd5_2-18-sone-dur-1300-ms';
    fname2suffix = 'JBS51-4544-Cd5_6-18-sone-dur-1300-ms';

    fname1   = [dir_where fname1suffix '.wav'];
    fname2   = [dir_where fname2suffix '.wav'];
end

legend1txt = 'signal1';
colour1 = 'b';

legend2txt = 'signal2';
colour2 = 'r';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Reading the wav files:
[signal1 fs1] = Wavread(fname1); 
[signal2 fs2] = Wavread(fname2); 

if fs1 ~= fs2
    error('Signals with different sampling frequency');
else
    fs = fs1;
end

[L idxL] = min([length(signal1) length(signal2)]);

hname = [];

t1 = (1:length(signal1))/fs;
t2 = (1:length(signal2))/fs;

rms1    = rmsdb(signal1)+100;
signal2 = setdbspl(signal2,rms1);

% noise  = noise1;
% fnamesuffix = fname1suffix;
ylimits = [-55 5];

if bPart1
    fc = 20;
    yenv1  = il_get_envelope(signal1,fs,fc);
    yenv2  = il_get_envelope(signal2,fs,fc);

    figure;
    Text_ylabel = 'Amplitude';
    Text_xlabel = 'Time [s]';

    subplot(2,1,1)
    plot(t1,signal1,colour1); hold on, grid on
    plot(t1,  yenv1,'k','LineWidth',2);
    plot(t1, -yenv1,'k','LineWidth',2);
    ylabel(Text_ylabel)
    xlabel(Text_xlabel)
    title(legend1txt)
    ha = gca;

    subplot(2,1,2)
    plot(t2, signal2,colour2); hold on, grid on
    plot(t2,   yenv2,'k','LineWidth',2);
    plot(t2,  -yenv2,'k','LineWidth',2);
    ylabel(Text_ylabel)
    xlabel(Text_xlabel)
    title(legend2txt)
    ha(end+1) = gca;

    plotOpts = [];
    plotOpts.I_Width = 11;
    plotOpts.FontSize = FontSize;

    h(end+1) = gcf;
    % hname{end+1} = ['waveform' pianoID];
    h_old = h(end); 
    h(end) = Figure2paperfigureT(h(end),2,1,plotOpts);
    close(h_old);

    figure;
    plot(t1,  yenv1, colour1,'LineWidth',2); hold on; grid on
    plot(t2,  yenv2, colour2,'LineWidth',2); 
    legend(legend1txt,legend2txt);
    ylabel(Text_ylabel)
    xlabel(Text_xlabel)
   
end

if bPart2

    fc = 20;
    yenv1  = il_get_envelope(signal1,fs,fc);
    yenv2  = il_get_envelope(signal2,fs,fc);
    
    plotOpts = [];
    plotOpts.I_Width = 10;
    plotOpts.FontSize = FontSize;
    plotOpts.I_KeepTicks = 0;

    figure;
    subplot(1,2,1)
    plot(t1,20*log10(abs(2*signal1)/2e-5),colour1); hold on
    ylabel('Sound Pressure Level [dB]');
    xlabel('Time [s]')
    title(legend1txt)

    ha = gca;

    plot(t1,20*log10((2*yenv1)/2e-5),'k','LineWidth',2);
    legend('Waveform','1. Hilb+LPF');

    %%%
    subplot(1,2,2)
    plot(t2,20*log10(abs(2*signal2)/2e-5),colour2); hold on
    plot(t2,20*log10((2*yenv2)/2e-5),'k','LineWidth',2);
    title(legend2txt)
    xlabel('Time [s]')

    ha(end+1) = gca;
    linkaxes(ha,'xy');
    ylim([40 90]);
    
    h(end+1) = gcf;
    h_old = h(end);
    h(end) = Figure2paperfigureT(h(end),1,2,plotOpts);
    close(h_old);
    
    figure;
    plot(t1, 20*log10((2*yenv1)/2e-5), colour1,'LineWidth',2); hold on; grid on
    plot(t2, 20*log10((2*yenv2)/2e-5), colour2,'LineWidth',2); 
    ylabel(Text_ylabel)
    xlabel(Text_xlabel)
end

if bPart3

    dir = [Get_TUe_paths('ICRA_Tobias')];
    addpath([dir delim 'Tools'])

    % [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,10e-3);
    % [max_0 idx_0] = max(RMS_se1);
    % Ni = find(t<t_se1(idx_0),1,'last');
    % Nf = Ni + round(500e-3*fs);

    % % Calculate LTAS
    [ltass1, fHz1] = calcLTAS(signal1,fs);
    [ltass2, fHz2] = calcLTAS(signal2,fs);

    dBMax = max( max([ltass1 ltass2]) );
    ltass1 = ltass1-dBMax;
    ltass2 = ltass2-dBMax;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    % subplot(1,2,1);
    semilogx(   fHz1,ltass1,colour1, ... 
                fHz1,ltass2,'-k'); hold on, grid on;
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
    % xlim([400 4200]);

    h(end+1) = gcf;
    plotOpts = [];
    plotOpts.I_Width = 11;
    plotOpts.FontSize = FontSize;
    h_old = h(end);
    h(end) = Figure2paperfigureT(h(end),1,1,plotOpts);
    % hname{end+1} = ['spectrum-' pianoID];
    close(h_old);

    legend({'piano','noise'},'location','northeast')

    rmpath([dir delim 'Tools'])  

end

disp('')
    
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

