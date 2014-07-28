function VoD_one_signal
% function VoD_one_signal
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 16/06/2014
% Last update: 17/06/2014 % Update this date manually
% Last used: 17/06/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDo_OB          = 0;
bDo_Gammatone   = 0;
bDo_STFT        = 1;

info.bSave      = 0;
info.wntype     = 1; % Hanning
% info.wntype     = 0; % Rectangular

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables:

h   = []; % empty handles
ha1 = []; % empty handles

nfft      = 4096*2; % N-point FFT and STFT
K         = nfft/2; % K effective FFT points
wlen_s    = 15e-3;  % STFT
overlap_p = 5;      % STFT

mode      = 2; % 2 to 5
fieldtype = 2; % 1 = far-field; 2 = near-field

paths.outputs   = Get_TUe_paths('outputs');
misc = Get_VoD_params(0,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading WAV-files and temporal plot

filemeasured = ['D:\MATLAB\Output\20140616-learning\modus-' num2str(mode-1) '_v2-' num2str(fieldtype) 'filt.wav']; % filt = LPF 
filemodelled = ['D:\MATLAB\Output\20140616-learning\modus-' num2str(mode-1) '-v_'  num2str(fieldtype) 'filt.wav'];
% filemodelled = ['D:\MATLAB\Output\20140616-learning\modus-' num2str(mode-1) '-v_'  num2str(fieldtype) '.wav'];
% filemodelled = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new\modus-1-v_2.wav';

% filename{1} = filemeasured;
% filename{2} = filemodelled;

[ynear fs]  = wavread(filemeasured);
ynearp      = wavread(filemodelled);

info.fs     = fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDo_STFT
    
    wlen = round(wlen_s * fs); % in samples
    
    figure;
    subplot(2,1,1)
    
    stft(ynear, fs, nfft, wlen, overlap_p, info.wntype);
    halocal = gca;
    
    cbYLim = [-120 -35];
    hcb = colorbar;
    set(hcb,'YLim',cbYLim);
    
    subplot(2,1,2)
    stft(ynearp, fs, nfft, wlen, overlap_p, info.wntype);
    halocal(end+1) = gca;
    hcb = colorbar;
    set(hcb,'YLim',cbYLim);
    
    linkaxes(halocal,'xy');
    
    ylim([0 1250])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = misc.t;
tp= ( 0:length(ynearp)-1 )/fs;

figure(1);
        
subplot(2,1,1)
plot(   t,ynear), grid on
legend('near-field') %,'far-field')
title(filemeasured)
haxis = gca;

subplot(2,1,2)
plot(   tp,ynearp), grid on
title(filemodelled)
haxis(end+1) = gca;

linkaxes(haxis,'xy');

samples = 20001;
ynear   = buffer(ynear,samples,0);
ynearp  = buffer(ynearp,samples,0);

ynear(:,end) = [];
ynearp(:,end) = [];
ynearp(:,1)=[];
ynear(:,1)=[];

RMS  = rmsdb(ynear);
RMSp = rmsdb(ynearp);

fprintf('Measured data, mean = %f; std = %f dB RMS\n',mean(RMS ),std(RMS ));
fprintf('Modelled data, mean = %f; std = %f dB RMS\n',mean(RMSp),std(RMSp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show_figures_one_by_one(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Windowing:
info.bNewFigure     = 0;
N = size(ynear,1);
M = size(ynear,2);
win                 = Get_window(info.wntype,N,M);
ynear               = ynear.*win;

N = size(ynearp,1);
M = size(ynearp,2);
[winp info.wtype]   = Get_window(info.wntype,N,M);
ynearp              = ynearp.*winp;

RMSw  = rmsdb(ynear); % RMS with window
RMSpw = rmsdb(ynearp);

fprintf('Measured data (windowed), mean = %f; std = %f dB RMS\n',mean(RMSw ),std(RMSw ));
fprintf('Modelled data (windowed), mean = %f; std = %f dB RMS\n',mean(RMSpw),std(RMSpw));

% OB analysis / one-third OB analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = [];
s = [];

lim_max = size(ynearp,2); % assuming size(Modelled) = size(Measured)
    
if bDo_OB
        
        figure(2);
                
        subplot(3,1,1)
        y2plot = [ynear ynearp];
        freqfft(y2plot,K,info);
         
        subplot(3,1,2)
        freqfft(ynear,K,info);
        title('measured')
        
        idx2plot = 1:lim_max;
        % idx2plot = 1:4;
        subplot(3,1,3)
        freqfft(ynearp(:,idx2plot),K,info);
        title('modelled')
        Print_text_on_plot(num2str( rmsdb(ynearp(:,idx2plot)) ),50);
         
        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);
         
        disp('Processing measured file');
        [Py_n1, Fc] = Filterbank_analysis( ynear,fs, 0);

        disp('Processing modelled file');
        [Ppy_n1, Fc] = Filterbank_analysis( ynearp,fs, 0);

        Delta = Ppy_n1 - Py_n1;
 
        m = [m, transpose( mean( transpose(Delta) ) )];
        s = [m, transpose(  std( transpose(Delta) ) )];

        [mean2plot,std2plot] = barweb_prepare_data(Py_n1,Ppy_n1);
 
        figure;
        Colors = [1 1 1;0.75 0.75 0.75];
        stPlot.Title = ['Field: near, Mode X: dB RMS per frequency band'];
        stPlot.XLabel = 'f_c [Hz]'; 
        stPlot.YLabel = 'Relative amplitude [dB]';
        stPlot.SeriesLabel = {'Measured','Modelled'};
 
        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        ha1(end+1) = gca;
        ha = gca;
        Set_XTick(ha,round(Fc));
        grid on

        ylabel('RMS value per band [dB]')
        xlabel('Central frequency f_c [Hz]')

        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('Processing measured file');
        [Pty_n1, Fc] = Filterbank_analysis( ynear,fs, 1);

        disp('Processing modelled file');
        Ptpy_n1 = Filterbank_analysis( ynearp,fs, 1);
 
        [mean2plot,std2plot] = barweb_prepare_data(Pty_n1,Ptpy_n1);
 
        figure;
        Colors = [1 1 1;0.75 0.75 0.75];
        stPlot.Title = ['Mode X: dB RMS per frequency band'];
        stPlot.XLabel = 'f_c [Hz]'; 
        stPlot.YLabel = 'Relative amplitude [dB]';
        stPlot.SeriesLabel = {'Measured','Modelled'};

        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        ha1(end+1) = gca;
        ha = gca;
        grid on

        Set_XTick(ha,round(Fc))

        ylabel('RMS value per band [dB]')
        xlabel('Central frequency f_c [Hz]')

        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);
        % Show_figures_one_by_one(0.5);
end

% Gammatone
if bDo_Gammatone

    disp('Processing measured file');
    [Pgy_n1 , Fc] = Filterbank_analysis( ynear,fs, 2);

    disp('Processing modelled file');
    [Pgpy_n1, Fc] = Filterbank_analysis( ynearp,fs, 2);

    [mean2plot,std2plot] = barweb_prepare_data(Pgy_n1,Pgpy_n1);

    figure;
    Colors = [1 1 1;0.75 0.75 0.75];
    stPlot.Title = ['Mode X: dB RMS per frequency band (Gammatone)'];
    stPlot.XLabel = 'f_c [Hz]'; 
    stPlot.YLabel = 'Relative amplitude [dB]';
    stPlot.SeriesLabel = {'Measured','Modelled'};

    barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
    ha1(end+1) = gca;
    ha = gca;
    Set_XTick(ha,round(Fc));
    grid on

    ylabel('RMS value per band [dB]')
    xlabel('Central frequency f_c [Hz]')
    
    h(end+1) = gcf;
    set(h(end),'Position', [1 41 1366 658]);

end

linkaxes(ha1,'y');
set( ha1(1),'YLim',[-80 -20]);

for j = 1:length(h)
    if info.bSave
        Saveas(h(j),[paths.outputs 'freq_analysis-' num2str(j)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end