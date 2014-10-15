function r20140930_FIA_KUL(options)
% function r20140930_FIA_KUL(options)
%
% 1. Description:
%       If assessment is already done, run this script with bAssess = 0
% 
%       Make sure all directories of NMT and NMTAddOns toolbox are added to
%       MATLAB path if you run this script.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 0;
%       r20140930_FIA_KUL(options);
%
% % 3.2 Example 2:
%       options.bSave = 1;
%       options.bAssess = 1; % Generates errors using stored CP810-simulated 
%                            % files (in quiet and in noise)
%       options.bDoLISTwhite = 1;
%       options.bDoLISTssn = 1;
%       r20140930_FIA_KUL(options);
% 
% % 3.3 Example 3:
%       options.bSave = 1;
%       options.bAssess = 0; % Results have to be previously generated
%       options.bDoLISTwhite = 1;
%       options.bDoLISTssn = 1;
%       options.bUseCleanSpeech = 1; % No CP810-simulation in F0ref
%       r20140930_FIA_KUL(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/09/2014
% Last update on: 07/10/2014 % Update this date manually
% Last use on   : 07/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

bDiary = 0;
Diary(mfilename,bDiary);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file in modules:
options = Ensure_field(options,'bDoPDA',1);     % Figure 4
options = Ensure_field(options,'bDoPlotSNR',1); % Figure 2, 3.a
options = Ensure_field(options,'bDoSpeech',1);  % Figure 5
options = Ensure_field(options,'bDoPlotF0',1); % Figure 3.b
options = Ensure_field(options,'bDoAdditionalAnalysis',0); % analysis as sent to TF on 07/10/2014
% Figure 1 was done in LibreOffice during KUL times

bDoPDA      = options.bDoPDA;
bDoPlotSNR  = options.bDoPlotSNR;
bDoSpeech   = options.bDoSpeech;
bDoPlotF0   = options.bDoPlotF0;
bDoAdditionalAnalysis = options.bDoAdditionalAnalysis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = Ensure_field(options,'bSave',1);
options = Ensure_field(options,'bDoLISTwhite',1);
options = Ensure_field(options,'bDoLISTssn',1);
options = Ensure_field(options,'bAssess',0);
options = Ensure_field(options,'bCleanSpeech',0);   % bCleanSpeech = 0: CP810 in quiet is reference
                                                    % bCleanSpeech = 1: wdz from LIST in quiet is reference
options = Ensure_field(options,'bDoNMT',0);
options = Ensure_field(options,'dest_folder_fig',[Get_TUe_paths('outputs') 'tmp-Figure-0417' delim]);
options = Ensure_field(options,'t_silence',0);
Mkdir(options.dest_folder_fig);

bDoLISTwhite    = options.bDoLISTwhite;
bDoLISTssn      = options.bDoLISTssn;
bAssess         = options.bAssess;

h = [];

%%
if bDoPDA == 1
    
    if bDoLISTwhite
        options.typenoise = 'white';
        options.bAssess = bAssess;
        [xx outWhite]   = experiment_report_20140228_Physical_validation_LIST(options);
    end

    if bDoLISTssn
        options.typenoise = 'SSN';
        options.bAssess = bAssess;
        [xx outSSN]     = experiment_report_20140228_Physical_validation_LIST(options);
    end

    close all

    %%
    % Figure 4.a
    figure
    stPlot.YLabel   = 'voiced error [%]';
    stPlot.yLim     = [0 100];
    stPlot.SeriesLabel  = {'white-noise','SSN'};
    Colors = [1 1 1; 0.6 0.6 0.6];
    h(end+1) = Plot_measure([outWhite.m1 outSSN.m1],[outWhite.s1 outSSN.s1], '', stPlot, Colors);
    xlim([0 5.5]); % to exclude SNR = -5 dB
    stPlot = [];

    % Figure 4.b
    figure
    stPlot.YLabel   = 'unvoiced error [%]';
    
    if options.bCleanSpeech == 0
        stPlot.yLim     = [0 50];
    else
        stPlot.yLim     = [0 70];
    end
    
    stPlot.SeriesLabel  = {'white-noise','SSN'};
    h(end+1) = Plot_measure([outWhite.m2 outSSN.m2],[outWhite.s2 outSSN.s2], '', stPlot, Colors);
    xlim([0 5.5]); % to exclude SNR = -5 dB
    stPlot = [];

    % Figure 4.c
    figure
    stPlot.YLabel    = 'gross error [%]';
    
    if options.bCleanSpeech == 0
        stPlot.yLim     = [0 20];
    else
        stPlot.yLim     = [0 25];
    end
    stPlot.SeriesLabel  = {'white-noise','SSN'};
    h(end+1) = Plot_measure([outWhite.m3 outSSN.m3],[outWhite.s3 outSSN.s3], '', stPlot, Colors);
    xlim([0 5.5]); % to exclude SNR = -5 dB
    stPlot = [];

    outSSN.t_total
    outSSN.t_unvoiced
    outSSN.t_voiced
    %%
    figure
    stPlot.YLabel   = '\Delta error [%]';
    stPlot.yLim     = [0 5];
    stPlot.SeriesLabel  = {'\Delta vErr','\Delta uvErr','\Delta gErr','\Delta Total Err'};
    Colors = [1 1 1; 0.7 0.7 0.7; 0.4 0.4 0.4; 0 0 0];
    h(end+1) = Plot_measure([outSSN.m11-outWhite.m11 outSSN.m12-outWhite.m12 outSSN.m13-outWhite.m13 outSSN.m99-outWhite.m99],[0*outWhite.s99 0*outWhite.s99 0*outWhite.s99 0*outWhite.s99], '', stPlot, Colors);
    xlim([0 5.5]); % to exclude SNR = -5 dB
    stPlot = [];


    disp(['Total time: ' num2str(outSSN.t_total) ' [s]'])
    disp('SSN error - White error / error * total time [%]')
    numDecimals = 2;
    tmp = [outSSN.m99 outWhite.m99 outSSN.m99-outWhite.m99 (outSSN.m99-outWhite.m99)*outSSN.t_total/100];
    tmp = Round(tmp,numDecimals);
    var2latex(tmp)

    if options.bSave
        for k = 1:length(h)
            if options.bCleanSpeech == 0
                Saveas(h(k),[options.dest_folder_fig 'LISTf-error-' num2str(k)]);
            else
                Saveas(h(k),[options.dest_folder_fig 'LISTf-error-' num2str(k) 'no-SP-simulation']);
            end
        end
    end
end

%% Plot SNR = 10 dB

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio-0417'  delim];
bGenerateFileSNR = 0; % change this value to 1 if some audio files are missing
tmp = Get_TUe_subpaths('db_speechmaterials');
root_folder = tmp.allfiles_LISTf;

if bDoPlotSNR
    bCreateWhitenoise = 0; % Whitenoise for LIST-f

    options.nAnalyser = 10;

    filename = [root_folder 'wdz2.wav']; % 'Elke zaterdag ga ik naar de markt'
    Wavread(filename);
    filenoise1 = [root_folder 'whitenoise-LISTf.wav'];
    filenoise2 = [root_folder 'wivineruis.wav'];

    options.CalMethod = 5;
    options.bPlot = 0;

    if bCreateWhitenoise
        [x1 fs1] = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
        x1RMS = rmsdb(x1);
        [x2 fs2] = Wavread([tmp.allfiles_PB 'whitenoise.wav']); % originally -10.78 dBFS RMS
        x2 = setdbspl(x2,x1RMS+100);
        rmsdb(x2)
        Wavwrite(x2,fs1,[tmp.allfiles_LISTf 'whitenoise-LISTf']);
    end

    out_file1 = PsySoundCL(filename,options);
    out_noise1 = PsySoundCL(filenoise1,options);
    out_noise2 = PsySoundCL(filenoise2,options);
    f = out_file1.f;

    SNR2plot = 10;
    
    figure; % Figure 2.1
    semilogx(   f, out_file1.DataSpecOneThirdAvg, 'bo--', ...
                f, out_noise1.DataSpecOneThirdAvg-SNR2plot, 'kx--', ...
                f, out_noise2.DataSpecOneThirdAvg-SNR2plot,'r>-');
    grid on, hold on
    % legend('wdz2, 65 dB SPL','white noise, 55 dB SPL','SSN, 55 dB SPL')
    legend('LIST-f: wdz2','white noise','SSN')
    xlabel('Frequency [Hz]')
    ylabel('Sound Pressure Level [dB]')
    ylim([20 65])
    xlim([50 22050])

    h(end+1) = gcf;

    %%%
    
    diff_noise1 = out_file1.DataSpecOneThirdAvg - (out_noise1.DataSpecOneThirdAvg-SNR2plot);
    diff_noise2 = out_file1.DataSpecOneThirdAvg - (out_noise2.DataSpecOneThirdAvg-SNR2plot);
    
    figure; % Figure 2.2
    semilogx(   f, diff_noise1, 'kx--', ...
                f, diff_noise2, 'r>-');
    legend('SNR using white noise','SNR using SSN')
    xlabel('Frequency [Hz]')
    ylabel('Estimated SNR [dB]')
    ylim([-20 30])
    xlim([50 22050])
    grid on, hold on
    semilogx([10 20000],[SNR2plot SNR2plot],'k--','LineWidth',2)
    
    h(end+1) = gcf;
    
    if options.bSave
        Saveas(h(end-1),[options.dest_folder_fig 'speech-SNR-' num2str(SNR2plot)]);
        Saveas(h(end  ),[options.dest_folder_fig 'speech-SNR-estimated']);
    end
    
    if bGenerateFileSNR
        [x1 fs1] = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
        x1RMS = rmsdb(x1);
        x1 = [x1; x1; x1];
        [x2 fs2] = Wavread(filename); % originally -10.78 dBFS RMS
        y = (From_dB(-SNR2plot)*x1(1:length(x2))) + x2;
        fileGenerated = [options.dest_folder 'wdz2-SNR-10-dB'];
        Wavwrite(y,fs1,fileGenerated);
    else
        fileGenerated = [options.dest_folder 'wdz2-SNR-10-dB.wav'];
        y = Wavread(fileGenerated);
    end

    if bGenerateFileSNR
        try
            f1 = ['CP810wdz2-10.wav'];
            [xx1 fss1] = Wavread([tmp.fda_eval_LISTf 'wav-CP810-SSN' delim f1]);

            yy1 = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);

            xx2 = resample(From_dB(5+rmsdb(yy1)-rmsdb(xx1))*xx1,44100,fss1);
            f1 = [options.dest_folder f1];
            Wavwrite(xx2,44100,f1);
        end
    else
        f1 = ['CP810wdz2-10.wav'];
        f1 = [options.dest_folder f1];
    end

    out_file1orig = PsySoundCL(fileGenerated,options);
    out_file1CP810 = PsySoundCL(f1 ,options);

    % Figure 3.a
    figure;
    semilogx(   f, out_file1orig.DataSpecOneThirdAvg, 'bo--'); hold on
    semilogx(   f, out_file1CP810.DataSpecOneThirdAvg,'kx-','LineWidth',2);
    grid on
    legend('wdz2','wdz2, CP810 sim.')
    xlabel('Frequency [Hz]')
    ylabel('Sound Pressure Level [dB]')

    ylim([20 65])
    xlim([50 9000])

    h(end+1) = gcf;

    if options.bSave
        Saveas(h(end),[options.dest_folder_fig 'speech-clean-CP810-' num2str(SNR2plot)]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%
if bDoPlotF0 
    outputfile  = [tmp.fda_eval_LISTf 'praat' delim 'wdz2.txt'];
    [tF0 F0]    = Get_F0_praat_from_txt(outputfile);

    outputfile  = [tmp.fda_eval_LISTf 'praat-CP810' delim 'CP810wdz2.txt'];
    [tF0test F0_CP810]    = Get_F0_praat_from_txt(outputfile);
    
    figure;
    plot(tF0, F0, '--'), hold on
    plot(tF0test-0.15, F0_CP810, 'k-','LineWidth',2);
    ylabel('F_0 [Hz]')
    xlabel('time [s]')
    grid on
    legend('wdz','wdz2, CP810 sim.')
    xlim([0 3.5])
    ylim([100 300])
    
    h(end+1) = gcf;
    if options.bSave
        Saveas(h(end),[options.dest_folder_fig 'speech-clean-CP810-F0-estimation']);
    end
end

%%
if bDoSpeech
    optionsSpeech = [];
    optionsSpeech.bDoSpeechTests = 1;
    optionsSpeech = Ensure_field(optionsSpeech,'bDoLT',1);
    optionsSpeech = Ensure_field(optionsSpeech,'bDoLL',0);
    optionsSpeech = Ensure_field(optionsSpeech,'bDoMT',0);

    % It is not working in R2013a...
    experiment_report_20131127_pooled_FIA2014(optionsSpeech);

    disp('')
end

if bDoAdditionalAnalysis
    r20141007_FIA_KUL_additional_analysis;
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = Plot_measure(m,s,title_str,stPlot, Colors)

if nargin < 4
    stPlot = [];
end

nTicks = 10;
stPlot.figPos       = [0 0 1024 300];
stPlot.xTickLabel   = {'Q','20','10','5','0','-5'};
stPlot.xTick        = 1:length(stPlot.xTickLabel);
stPlot.xLim         = [ 0   max(stPlot.xTick)+1];
stPlot              = Ensure_field(stPlot, 'yLim', [0 50]);
stepTicks           = ceil( (stPlot.yLim(2)-stPlot.yLim(1))/nTicks );

stPlot.yTick        = [stPlot.yLim(1)+stepTicks:stepTicks:stPlot.yLim(2)-stepTicks];
stPlot.Title        = title_str;
stPlot.XLabel       = 'SNR (dB)';
stPlot              = Ensure_field(stPlot, 'YLabel','%');

if size(m,2) == 1
    
    if ~exist('Colors','var')
        Colors              = [1 1 1];  
    end
    m = [m 0*m]; % trick to plot bars series with x-axis in steps of 1
    s = [s 0*s];
else % then size == 2
    stPlot  = Ensure_field(stPlot,'SeriesLabel',{'PB male','PB female','LIST-f'});
    if ~exist('Colors','var')
        Colors              = [1 1 1; 0.75 0.75 0.75; 0 0 0];  
    end
end

if isfield(stPlot, 'SeriesLabel')
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
else
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors);
end

grid on, hold on

h = gcf;
set(h,'Position', stPlot.figPos);

Handle = gca;
set(Handle,'YTick',stPlot.yTick)
set(Handle,'XTick',stPlot.xTick)
set(Handle,'XLim',stPlot.xLim)
set(Handle,'YLim',stPlot.yLim)
set(Handle,'XTickLabel',stPlot.xTickLabel)

set(h, 'PaperPositionMode','auto')