function Run_HRIR
% function Run_HRIR
%
% Uses Oldenburg's database of HRIR measured in 4 channels per ear
% MAT-files:
%   Nomenclature:   anechoic_distcm_300_el_0_az_-45.mat
%   It contains:
%       - azimuth
%       - data (N x 8 channels):    Ch 1 and 2 - in-ear Left and Right
%                                   Ch 3 and 4 - BTE front Left and Right
%                                   Ch 5 and 6 - BTE mid Left and Right
%                                   Ch 7 and 8 - BTE rear Left and Right
%                            Note:  N is 4800, then audio files are 100 ms long
%       - distance_cm:              300 cm = 3  m = far-field
%                                    80 cm = 0.8m = near-field
%       - elevation
%       - environment ('Anechoic')
%       - fs
%       - nBits
%       - nChannels
%       - nSamples
% 
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 08/05/2014
% Last update: 15/05/2014 % Update this date manually
% Last used: 15/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
info.bSave = 0;

h = []; % we initialise handle for Figures
paths = Get_TUe_paths;

% We are going to work with one HRIR: far field, anechoic room, elevation of 
%   0 degrees and azimuth of -45 degrees:
ambient     = 'anechoic';
elevation   = 0;
Az          = [0 -45];
dist        = 300; % 300 cm = far field
dB_SPL      = 65; % reference: Left audio file

for count = 1:length(Az)
    azimuth     = Az(count);
    
    % Then we will load the stored data �
    %   (Nomenclature: anechoic_distcm_300_el_0_az_-45.mat)
    filename = [ambient '_distcm_' num2str(dist) '_el_' num2str(elevation) '_az_' num2str(azimuth) '.mat'];
    load([paths.db_HRIR_Oldenburg 'HRIR_database_mat' delim 'hrir' delim ambient delim filename]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot of frequency responses
    K_factor    = From_dB(+35.2); % Visual assessment to match freq response to paper 
    in_ear      = K_factor*data(:,1:2); % Only in-ear channels

    N = 4096*2; % N-FFT points
    K = N/2;

    info.fs = fs;
    freqfft(in_ear,K,info); % just to generate plot
    xlim([100 11000])
    ylim([-30 30])
    legend('in-ear L','in-ear R')

    [H HdB f] = freqfft(in_ear,K,info);
    
    % % Smoothed FFT spectrum, uncomment following 6 lines to perform it
    % smoothfft(HdB, f); % just to generate (smoothed) plot
    % xlim([100 11000])
    % ylim([-30 30])
    % legend('in-ear L','in-ear R')
    % HdB_smooth = smoothfft(HdB, f);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convolution:
    t_duration = 100e-3;
    Naudio = round(info.fs*t_duration);
    % White noise:
    filename = [paths.outputs 'gwn_test1'];
    if info.bSave
        y   = wgn(Naudio,1,1);
        y   = setdbspl(y,dB_SPL);

        wavwrite(y,fs,filename);
    else
        [y FFs] = wavread(filename);
    end

    t = ( 1:length(y) )/info.fs;
    t = t(:);

    [yconvL tconv]  = Conv(in_ear(:,1), y, info);
    [yconvR]        = Conv(in_ear(:,2), y, info);

    temporal        = [in_ear(:,1) y];
    temporalConv    = [yconvL yconvR];
    freqfft(temporal,K,info);
    legend('in-ear L','white'); 

    h(end+1) = gcf;

    freqfft(temporalConv,K,info);
    legend('white-L-ear','white-R-ear')

    h(end+1) = gcf;

    % Lindemann's model
    sig = [yconvL yconvR]; % Binaural signal

    % Calculate binaural cross-correlation
    [cc,t,~,~,params] = lindemann1986(sig,fs,'T_int',6);

    % Set title string for the plot
    tstr = sprintf(['Binaural white noise, (az = %i, fs = %i Hz)\n' params.text],azimuth,fs);
    % Plot frequency channel 11, due to round(freqtoerb(500))==11
    f2plot = 1000;
    figure;
    plotlindemann1986(cc,t,'fc',f2plot,'title',name2figname(tstr));

    h(end+1) = gcf;
    set(h(end),'Position', [680 276 926 702])
    
    if info.bSave
        Wavwrite(yconvL,fs,[filename 'az' num2str(azimuth) '-L']);
        Wavwrite(yconvR,fs,[filename 'az' num2str(azimuth) '-R']);
    end
    
    % 4. Dau 1997
    insig = yconvL;
    [outsig, fc,~,allouts] = dau1997preproc(insig,fs);
    i = 10; % band number 10
    
    figure
    subplot(2,1,1)
    plot(   1:length(allouts.out01_filterbank(:,i)), allouts.out01_filterbank(:,i) ,...
            1:length(allouts.out02_ihc(:,i))       , allouts.out02_ihc(:,i))
	legend('gammatone filter', 'half-wave rectification')
    xlabel('[Sample number]')
    title(sprintf('Dau''s model for gammatone filter, fc = %.0f Hz',fc(i)))
    grid on
    
    subplot(2,1,2)
    plot(   1:length(allouts.out03_adaptloop(:,i)) , allouts.out03_adaptloop(:,i))
	legend('adaptive loops')
    xlabel('[Sample number]')
    grid on
    
end

if info.bSave
    for i=1:length(h)
        Saveas(h(i), [filename '-handle-' num2str(i)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename])
