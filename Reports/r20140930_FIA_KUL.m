function r20140930_FIA_KUL(options)
% function r20140930_FIA_KUL(options)
%
% 1. Description:
%       If assessment is done, run this script with bAssess = 0
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/09/2014
% Last update on: 03/09/2014 % Update this date manually
% Last use on   : 03/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

bDiary = 0;
Diary(mfilename,bDiary);

options = Ensure_field(options,'bSave',1);
bDoLISTwhite    = 1;
bDoLISTssn      = 1;
bAssess         = 0;

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

h = [];
close all

figure
stPlot.YLabel   = 'voiced error [%]';
stPlot.yLim     = [0 100];
stPlot.SeriesLabel  = {'white-noise','SSN'};
Colors = [1 1 1; 0.6 0.6 0.6];
h(end+1) = Plot_measure([outWhite.m1 outSSN.m1],[outWhite.s1 outSSN.s1], '', stPlot, Colors);
xlim([0 5.5]); % to exclude SNR = -5 dB
stPlot = [];

figure
stPlot.YLabel   = 'unvoiced error [%]';
stPlot.yLim     = [0 50];
stPlot.SeriesLabel  = {'white-noise','SSN'};
h(end+1) = Plot_measure([outWhite.m2 outSSN.m2],[outWhite.s2 outSSN.s2], '', stPlot, Colors);
xlim([0 5.5]); % to exclude SNR = -5 dB
stPlot = [];

figure
stPlot.YLabel    = 'gross error [%]';
stPlot.yLim     = [0 20];
stPlot.SeriesLabel  = {'white-noise','SSN'};
h(end+1) = Plot_measure([outWhite.m3 outSSN.m3],[outWhite.s3 outSSN.s3], '', stPlot, Colors);
xlim([0 5.5]); % to exclude SNR = -5 dB
stPlot = [];

% figure
% stPlot.YLabel    = 'f0Dev [\Delta Hz]';
% stPlot.SeriesLabel  = {'white-noise','SSN'};
% h(end+1) = Plot_measure([outWhite.m4 outSSN.m4],[outWhite.s4 outSSN.s4], '', stPlot);
% xlim([0 5.5]); % to exclude SNR = -5 dB
% stPlot = [];

outSSN.t_total = outSSN.t_total-350*0.5;

figure
stPlot.YLabel   = '\Delta error [%]';
stPlot.yLim     = [0 5];
stPlot.SeriesLabel  = {'\Delta vErr','\Delta uvErr','\Delta gErr','\Delta vErr+uvErr+gErr'};
Colors = [1 1 1; 0.7 0.7 0.7; 0.4 0.4 0.4; 0 0 0];
h(end+1) = Plot_measure([outSSN.m11-outWhite.m11 outSSN.m12-outWhite.m12 outSSN.m13-outWhite.m13 outSSN.m99-outWhite.m99],[0*outWhite.s99 0*outWhite.s99 0*outWhite.s99 0*outWhite.s99], '', stPlot, Colors);
xlim([0 5.5]); % to exclude SNR = -5 dB
stPlot = [];


disp(['Total time: ' num2str(outSSN.t_total) ' [s]'])
disp('SSN error - White error / error * total time [%]')
disp(num2str([outSSN.m99 outWhite.m99 outSSN.m99-outWhite.m99 (outSSN.m99-outWhite.m99)*outSSN.t_total/100]))

if options.bSave
    for k = 1:length(h)
        Saveas(h(k),[Get_TUe_paths('outputs') 'LISTf-error-' num2str(k)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot SNR = 10 dB

bCreateWhitenoise = 0; % Whitenoise for LIST-f

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
options.nAnalyser = 10;

tmp = Get_TUe_subpaths('db_speechmaterials');
root_folder = tmp.allfiles_LISTf;
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
figure;
semilogx(   f, out_file1.DataSpecOneThirdAvg, 'bo--', ...
            f, out_noise1.DataSpecOneThirdAvg-SNR2plot, 'kx--', ...
            f, out_noise2.DataSpecOneThirdAvg-SNR2plot,'r>-');
grid on
legend('wdz2, 65 dB SPL','white noise, 55 dB SPL','SSN, 55 dB SPL')
xlabel('Frequency [Hz]')
ylabel('Sound Pressure Level [dB]')
ylim([20 65])
xlim([50 22050])

h(end+1) = gcf;

if options.bSave
    Saveas(h(end),[Get_TUe_paths('outputs') 'speech-SNR-' num2str(SNR2plot)]);
end

bGenerateFileSNR = 0;
if bGenerateFileSNR
    [x1 fs1] = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
    x1RMS = rmsdb(x1);
    x1 = [x1; x1; x1];
    [x2 fs2] = Wavread(filename); % originally -10.78 dBFS RMS
    y = (From_dB(-SNR2plot)*x1(1:length(x2))) + x2;
    fileGenerated = [Get_TUe_paths('outputs') 'wdz2-SNR-10-dB'];
    Wavwrite(y,fs1,fileGenerated);
else
    fileGenerated = [Get_TUe_paths('outputs') 'wdz2-SNR-10-dB.wav'];
    y = Wavread(fileGenerated);
end

if bGenerateFileSNR
    try
        f1 = ['CP810wdz2-10.wav'];
        [xx1 fss1] = Wavread([tmp.fda_eval_LISTf 'wav-CP810-SSN' delim f1]);
        
        yy1 = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
                
        xx2 = resample(From_dB(5+rmsdb(yy1)-rmsdb(xx1))*xx1,44100,fss1);
        f1 = [Get_TUe_paths('outputs') f1];
        Wavwrite(xx2,44100,f1);
    end
else
    f1 = ['CP810wdz2-10.wav'];
    f1 = [Get_TUe_paths('outputs') f1];
end

out_file1orig = PsySoundCL(fileGenerated,options);
out_file1CP810 = PsySoundCL(f1 ,options);

figure;
semilogx(   f, out_file1orig.DataSpecOneThirdAvg, 'bo--', ...
            f, out_file1CP810.DataSpecOneThirdAvg, 'kx-' )
grid on
legend('wdz2, LIST-f','wdz2, CP810 sim.')
xlabel('Frequency [Hz]')
ylabel('Sound Pressure Level [dB]')

ylim([20 65])
xlim([50 9000])

h(end+1) = gcf;

if options.bSave
    Saveas(h(end),[Get_TUe_paths('outputs') 'speech-clean-CP810-' num2str(SNR2plot)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
optionsSpeech = [];
optionsSpeech.bDoSpeechTests = 1;
optionsSpeech = Ensure_field(optionsSpeech,'bDoLT',1);
optionsSpeech = Ensure_field(optionsSpeech,'bDoLL',0);
optionsSpeech = Ensure_field(optionsSpeech,'bDoMT',0);

% It is not working in R2013a...
experiment_report_20131127_pooled_FIA2014(optionsSpeech);

disp('')
 
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