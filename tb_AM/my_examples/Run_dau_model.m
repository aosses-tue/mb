function Run_dau_model
% function Run_dau_model
%
% 1. Description:
% 
% 2. Additional info:
% 
% 3. Stand-alone example:
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/05/2014
% Last update on: 28/05/2014 % Update this date manually
% Last use on   : 07/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
info.bSave      = 1;
info.bGenerate  = 0;

signaltype = 0; % 0 = white noise; 1 = sine at 600 Hz

h = []; % we initialise handle for Figures
paths.outputs   = Get_TUe_paths('outputs');

dB_SPL          = 65; % reference: Left audio file

N   = 4096*2; % N-FFT points
K   = N/2;
fs  = 48000;

info.fs = fs;
info.typeplot = 2; % Linear scaled

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution:
t_silence_before= 100e-3;
t_duration      = 200e-3;
t_silence_after = 300e-3;

Nsil_bef    = round(info.fs*t_silence_before);
Nnoise      = round(info.fs*t_duration);
Nsil_aft    = round(info.fs*t_silence_after);

switch signaltype
    case 0
        filename = [paths.outputs 'dau_gwn_test1'];
        title1 = 'White noise';
        ymin = -0.15;
        ymax =  0.15;
        yminMU = -100;
        ymaxMU = 1500;
    case 1
        filename = [paths.outputs 'sine_600Hz'];
        title1 = 'Sine wave';
        ymin = -1;
        ymax =  1;
        yminMU = -100;
        ymaxMU = 1500;
end

if info.bGenerate
    
    switch signaltype
        case 0
            y = wgn(Nnoise,1,1);
        case 1
            y = Create_sin(600,t_duration,fs);
    end
    y   =  y(:);
    y   = [zeros(Nsil_bef,1); y; zeros(Nsil_aft,1)]; % silence at the beginning and at the end
    y   = setdbspl(y,dB_SPL);

    if info.bSave
        Wavwrite(y,fs,filename);
    end
else
    try
        [y fs] = wavread(filename);
        % info.fs has to be equal to fs
    catch
        warning('Set info.bGenerate to 1 and info.bSave to 1 and re-run this script. Maybe this will solve the problem...')
        error(['Wav file not found make sure it is in: ' paths.output])
    end
end

t = ( 1:length(y) )/info.fs;
t = t(:);

plot(t,y); grid on
Add_pdf2plot(y,gcf);
xlabel('Time [s]')
ylabel('Amplitude')
ylim([ymin ymax])
h(end+1)=gcf;

freqfft(y,K,info);
legend(title1)
xlim([0 8000])
ylim([-30 20])
xlabel('Time [s]')
ylabel('Amplitude')
h(end+1)=gcf;

% 4. Dau 1997
insig = y;

[outsig, fc,xx,allouts] = dau1997preproc(insig,fs);

haxis = [];

for i = 1:10:31 % band number 10 and 20
    figure
    subplot(2,1,1)
    plot(   t, allouts.out01_filterbank(:,i) ,...
            t, allouts.out02_ihc(:,i))
    legend('gammatone filter', 'half-wave rectification')
    xlabel('Time [s]')
    ylim([ymin ymax]/10)
    title(sprintf('Dau''s model for gammatone filter, fc = %.0f Hz',fc(i)))
    grid on
    haxis(end+1) = gca;
    
    subplot(2,1,2)
    plot(   t , allouts.out03_adaptloop(:,i))
    legend('adaptive loops')
    xlabel('Time [s]')
    grid on

    h(end+1)=gcf;
    haxis(end+1) = gca;
end
linkaxes(haxis,'x');
xlim([0 0.5])

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

if info.bSave
    for i=1:length(h)
        Saveas(h(i), [filename '-handle-' num2str(i)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename])