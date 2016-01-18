function r20160125_update_BATWOMAN
% function r20160125_update_BATWOMAN
%
% 1. Description:
%       The ICRA noises were generated using a previous version of this
%       script (method in icra5_noise4piano.m).
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Original file name: r20151223_update_Antoine.m
% Created on    : 11/01/2016
% Last update on: 11/01/2016 
% Last use on   : 11/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

bSave   = 1;
hFig    = [];
bCreate = 0;
nSig2plot = 2;

dir_out   = [Get_TUe_data_paths('lx_Text') 'lx2015-12-21-update-ICRA-Antoine' delim 'piano-sounds-new' delim];
Mkdir(dir_out);

dir = [Get_TUe_paths('ICRA_Tobias')];
addpath([dir delim 'Tools']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sounds to be processed:
note_test = 'A4';
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];

fname0suffix = 'NS19-A4_2';
fname1suffix = 'GRAF28-A4_3';
fname2suffix = 'JBS51-4544-A4_4';
fname0   = [dir_where fname0suffix '.wav'];
fname1   = [dir_where fname1suffix '.wav'];
fname2   = [dir_where fname2suffix '.wav'];

legend0txt = 'NS19 (A4)';
legend0txtn = 'ICRA noise from NS19 (A4)';
colour0 = 'b';

legend1txt = 'Graf28 (A4)';
legend1txtn = 'ICRA noise from Graf28';

legend2txt = 'JBS51-4544 (A4)';
legend2txtn = 'ICRA noise from JBS51-4544';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Reading the wav files:
[signal0 fs] = Wavread(fname0); 
[signal1 fs] = Wavread(fname1); 
[signal2 fs] = Wavread(fname2); 

t0 = ( 1:length(signal0) )/fs;
t1 = ( 1:length(signal1) )/fs;
t2 = ( 1:length(signal2) )/fs;

% 2. Aligning the waveforms:
[RMS_se0 t_se0] = rmsdb_sec(signal0,fs,10e-3);
[RMS_se1 t_se1] = rmsdb_sec(signal1,fs,10e-3);
[RMS_se2 t_se2] = rmsdb_sec(signal2,fs,10e-3);

[max_0 idx_0] = max(RMS_se0);
[max_1 idx_1] = max(RMS_se1);
[max_2 idx_2] = max(RMS_se2);

signal0 = From_dB(-10-8) * signal0;
signal1 = From_dB(-10  ) * signal1;
signal2 = From_dB(-10-8) * signal2;

idx_0 = find(t0 <= t_se0(idx_0),1,'last');
idx_1 = find(t1 <= t_se1(idx_1),1,'last');
idx_2 = find(t2 <= t_se2(idx_2),1,'last');

samples_diff = abs(idx_2 - idx_1);
if idx_1 < idx_2 
% then idx_2 is after
   [signal2 t2] = Do_alignment(t2,signal2,t2(samples_diff));
else
% then idx_2 is after
   [signal1 t1] = Do_alignment(t1,signal1,t1(samples_diff));
end

samples_diff = abs(idx_0 - idx_1);
if idx_1 < idx_0
   [signal0 t0] = Do_alignment(t0,signal0,t0(samples_diff)); 
else
    error('Waveforms aligned only assuming signal1 as reference')
end
    
[L idxL] = min([length(signal0) length(signal1) length(signal2)]);
signal0 = Do_truncate(signal0,L);
signal1 = Do_truncate(signal1,L);
signal2 = Do_truncate(signal2,L);

dur_ramp = 150; % ms
if idxL == 1
% then signal1 and signal2 were truncated
    signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp);
    signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp);
elseif idxL == 2
% then signal0 and signal2 were truncated
    signal0 = Do_cos_ramp(signal0,fs,0,dur_ramp);
    signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp);
elseif idxL == 3
    % then signal0 and signal2 were truncated
    signal0 = Do_cos_ramp(signal0,fs,0,dur_ramp);
    signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp);
end

clear t0 t1 t2;
t = (1:L) / fs;

% 2.2. Storing the aligned waveforms:
if bSave
    fname_out = sprintf('%s%s-att.wav',dir_out,fname0suffix);
    Wavwrite(signal0,fs,fname_out);
    fname_out = sprintf('%s%s-att.wav',dir_out,fname1suffix);
    Wavwrite(signal1,fs,fname_out);
    fname_out = sprintf('%s%s-att.wav',dir_out,fname2suffix);
    Wavwrite(signal2,fs,fname_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method    = 3; % i = 3 = ERB;
fname_out0 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname0suffix,method);
fname_out1 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname1suffix,method);
fname_out2 = sprintf('%snoise-%s-ICRA-meth-%.0f.wav',dir_out,fname2suffix,method);

if bCreate
    noise0   = icra5_noise4piano(signal0,fs,method);
    Wavwrite(noise0,fs,fname_out0);
    noise1   = icra5_noise4piano(signal1,fs,method);
    Wavwrite(noise1,fs,fname_out1);
    noise2   = icra5_noise4piano(signal2,fs,method);
    Wavwrite(noise2,fs,fname_out2);
else
    noise0  = Wavread(fname_out0);
    noise1  = Wavread(fname_out1);
    noise2  = Wavread(fname_out2);
end

switch nSig2plot
    case 0
        signal = signal0;
        noise  = noise0;
        fnamesuffix = fname0suffix;
        ylimits = [-75 -20];
        legendtxt = legend0txt;
        legendtxtn = legend0txtn;
    case 1
        signal = signal1;
        noise  = noise1;
        fnamesuffix = fname1suffix;
        ylimits = [-80 -25];
        legendtxt = legend1txt;
        legendtxtn = legend1txtn;
    case 2
        signal = signal2;
        noise  = noise2;
        fnamesuffix = fname2suffix;
        ylimits = [-75 -20];
        legendtxt = legend2txt;
        legendtxtn = legend2txtn;
end
fc = 20;
yenv  = il_get_envelope(signal,fs,fc);
yenvn = il_get_envelope(noise,fs,fc);

figure;
Text_ylabel = 'Amplitude';
Text_xlabel = 'Time [s]';

subplot(2,1,1)
plot(t,signal,colour0); hold on, grid on
plot(t,  yenv,'k','LineWidth',2);
plot(t, -yenv,'k','LineWidth',2);
ylabel(Text_ylabel)
title(legendtxt)
ha = gca;
 
subplot(2,1,2)
plot(t, noise,colour0); hold on, grid on
plot(t,  yenvn,'k','LineWidth',2);
plot(t, -yenvn,'k','LineWidth',2);
ylabel(Text_ylabel)
title(legendtxtn)
ha(end+1) = gca;

linkaxes(ha,'xy');
xlabel(Text_xlabel)
ylim([-0.049 0.049])
xlim([0 max(t)])

hFig(end+1) = gcf;
hFig(end) = Figure2paperfigure(hFig(end),2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate LTAS
[ltassP, fHz] = calcLTAS(signal,fs);
[ltassN     ] = calcLTAS(noise ,fs);

% Calculate MTFs
[mtfP,cfModHz] = calcMTF(signal,fs);
[mtfN        ] = calcMTF(noise ,fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% subplot(1,2,1);
semilogx(   fHz,ltassP,colour0, ... 
            fHz,ltassN,'-k'); hold on, grid on;
xlabel('Frequency [Hz]')
ylabel('LTAS [dB]')
legend({'signal','noise'},'location','southwest')
axis tight;
XTick = [50 100 250 500 1000 2000 4000 8000];
set(gca,'XTick'     ,XTick);
set(gca,'XTickLabel',XTick);
ha0 = gca;    
title(name2figname(fnamesuffix));

linkaxes(ha0,'xy');
ylim(ylimits)
xlim([80 8500]);

hFig(end+1) = gcf;
hFig(end) = Figure2paperfigure(hFig(end),1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise01 = 0.5*noise0+0.5*noise1;
noise02 = 0.5*noise0+0.5*noise2;
noise12 = 0.5*noise1+0.5*noise2;

if bSave
    fname_out01 = sprintf('%snoise-mixed-ICRA-meth-%.0f-01.wav',dir_out,method);
    Wavwrite(noise01,fs,fname_out01);
    fname_out02 = sprintf('%snoise-mixed-ICRA-meth-%.0f-01.wav',dir_out,method);
    Wavwrite(noise01,fs,fname_out02);
    fname_out12 = sprintf('%snoise-mixed-ICRA-meth-%.0f-01.wav',dir_out,method);
    Wavwrite(noise01,fs,fname_out12);
end

if bSave
    for i = 1:length(hFig)
        fname = sprintf('%sICRA-%.0f-fig-%.0f',dir_out,nSig2plot,i);
        Save_figure_as(hFig(i),fname,'epsc');
        Save_figure_as(hFig(i),fname,'emf');
        Save_figure_as(hFig(i),fname,'fig');
    end
end

rmpath([dir delim 'Tools'])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function yenv = il_get_envelope(insig,fs,fc)

if nargin < 3
    fc = 20;
end

yin = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv = filtfilt(b,a,yin);
