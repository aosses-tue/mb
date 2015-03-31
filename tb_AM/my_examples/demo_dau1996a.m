function demo_dau1996a(options)
% function demo_dau1996a(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       demo_dau1996a;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 07/10/2014 % Update this date manually
% Last use on   : 07/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    options = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot save/options:
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal options:
options = Ensure_field(options, 'dB_SPL', 77);
options = Ensure_field(options, 'dB_SPL_test', 85);

infilename         = ['dau1996a_noisemasker-' num2str(options.dB_SPL)     ];
infilename2        = ['dau1996a_testsignal-'  num2str(options.dB_SPL_test)];

h = []; % we initialise handle for Figures
paths.outputs   = Get_TUe_paths('outputs');

try
    [insig   fs] = Wavread([paths.outputs infilename  '.wav']);
    [insig2  fs] = Wavread([paths.outputs infilename2 '.wav']);
    options.bGenerate   = 0;
catch
    options.bGenerate   = 1;
    fs                  = 48000;
end

dB_SPL  = options.dB_SPL; % reference: Left audio file
dB_SPL2 = options.dB_SPL_test;

N   = 4096*2; % N-FFT points
K   = N/2;

options.fs = fs;
options.typeplot = 2; % Linear scaled

t_silence_before= 100e-3;
t_duration      = 200e-3;
t_silence_after = 300e-3;
t_total_duration = t_silence_before + t_duration + t_silence_after;

Nsil_bef    = round(options.fs*t_silence_before);
Nnoise      = round(options.fs*t_duration);
Nsil_aft    = round(options.fs*t_silence_after);

%% Generating the signal
filename    = [paths.outputs infilename];
filename2   = [paths.outputs infilename2];
title1 = 'White noise';
ymin = -0.15;
ymax =  0.15;
yminMU = -100;
ymaxMU = 1500;
    
if options.bGenerate
    
    % Gen1: white noise, band-pass filtered
    y = wgn(Nnoise,1,1);
    
    y   =  y(:); % ensures it is a column vector
    y   = setdbspl(y,dB_SPL);
    y   = [zeros(Nsil_bef,1); y; zeros(Nsil_aft,1)]; % silence at the beginning and at the end
    
    Wn = [20 5000]/(options.fs/2); % Normalised cutoff frequency        
    [b,a] = butter(4,Wn); % 8th-order
    insig = filtfilt(b,a,y); % Linear-phase implementation
        
    Wavwrite(insig,fs,filename);
    
    % Gen2: sine tone
    
    f       = 3000;
    dur     = 10e-3;
    win     = 1; % 1 = Hanning window
    [y2, t2]= Create_sin(f,dur,fs,win);
    y2      = setdbspl(y2,dB_SPL2);
    insig2      = [Gen_silence(0.2,fs); y2; Gen_silence(t_total_duration-max(t2)-0.2-1/fs,fs)];

    Wavwrite(insig2,fs,filename2);
    
end

t = ( 1:length(insig) )/options.fs;
t = t(:);

% Dau1996compare(insig,insig2,fs);

%% Similar to Dau1996a, Figure 4
figure;
plot(t,insig); grid on
Add_pdf2plot(insig,gcf);
xlabel('Time [s]')
ylabel('Amplitude')
ylim([ymin ymax])
h(end+1)=gcf;

%% Similar to Dau1996a, Figure 5
freqfft(insig,K,options);
legend(title1)
xlim([0 8000])
ylim([-30 20])
xlabel('Time [s]')
ylabel('Amplitude')
h(end+1)=gcf;

ti = 100e-3;
tf = 300e-3;
stats1 = Get_stats_from_audio(insig       ,ti,tf,fs);
stats2 = Get_stats_from_audio(insig+insig2,ti,tf,fs);
stats3 = Get_stats_from_audio(insig2      ,ti,tf,fs);

figure;
plot(   stats1.pdf_centres,stats1.pdf, ...
        stats2.pdf_centres,stats2.pdf);
legend('noise alone','noise + signal')
%% Processing 

[outsig , fc ,allouts ] = dau1996preproc(insig         ,fs);
[outsig2, fc2,allouts2] = dau1996preproc(insig + insig2,fs);

idx     = max(find(fc<3000));
fcentre = fc(idx);

haxis = [];

%% Similar to Dau1996a, Figure 6
% Apparently this filter bank is the replaced one (see Dau1997b)
figure;
freqfft(allouts.out01_filterbank(:,idx) ,K,options);
legend(title1)
xlim([0 8000])
ylim([-30 20])
xlabel('Time [s]')
ylabel('Amplitude')
h(end+1)=gcf;

%% Similar to Dau1996a, Figure 7
for i = idx 
    
    figure
    subplot(3,1,1)
    plot(   t, allouts.out04_LPF(:,idx) ); hold on
    
    subplot(3,1,2)
    plot(   t, allouts2.out04_LPF(:,idx) ); hold on
    
    D = allouts2.out04_LPF(:,idx)-allouts.out04_LPF(:,idx);
    subplot(3,1,3)
    plot(   t,  25e-4*D); hold on
    
    h(end+1)=gcf;
    haxis = gca;
end
linkaxes(haxis,'x');
% xlim([0 0.5])

figure; 
mesh(allouts.out03_adaptloop)
xlabel('Band number'); 
ylabel('Sample number'); 
zlabel('Amplitude [MU]'); 
h(end+1)=gcf;
ha = gca;

zlim([yminMU ymaxMU])
set(ha,'CameraPosition',[79.9938 138832 22702.7]);
colorbar('vert')
% 	CameraPosition = [79.9938 138832 22702.7]

if options.bSave
    for i=1:length(h)
        Saveas(h(i), [filename '-handle-' num2str(i)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
