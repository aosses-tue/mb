function experiment_report_20140402_compare_F0ref(info)
% function experiment_report_20140402_compare_F0ref(info)
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
% Last update on: 09/10/2014
% Last use on   : 09/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    info = [];
end

info = Ensure_field(info,'bPB_male'  ,1);
info = Ensure_field(info,'bPB_female',1);
info = Ensure_field(info,'bLISTf_SSN',1);

bPB_male = info.bPB_male;
bPB_female = info.bPB_female;
% bLISTf_white = 0;
bLISTf_SSN = info.bLISTf_SSN;

info = Ensure_field(info,'bAssess',0);
info.bAnalyse   = 1;
info.isPaper    = 1;
info = Ensure_field(info,'bSave', 0);

if info.bSave
    disp([mfilename '.m: Figures and other files will be stored and might replace other files. Press any button to continue'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause()
else
    disp([mfilename '.m: Figures and other files won''t be stored...change bSave to 1 if want to save results to file...'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause(3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info_PB = info;
info_PB.F0reference     = 'laryngeal';
info_PB.F0max           = 400;

info_PB.bCleanSpeech    = 1; % if 0: CP810 simulations will be used
info_PB.bPlot           = 0;
info_PB.speaker         = 'sb'; % Female 
info_PB.results_dir_name = 'results';
info_PB.figures_dir_name = 'figures';


Colors = [0 0 0; 1 1 1];
h = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PB database
if bPB_female
    [info_PB tmpPBf] = experiment_report_20140402_Physical_validation(info_PB);
    close all
end

if bPB_male
    info_PB_male = info_PB;
    info_PB_male.speaker = 'rl'; % Male

    [info_PB_male tmpPBm] = experiment_report_20140402_Physical_validation(info_PB_male);
    close all
end

if bPB_male & bPB_female
    % Analysis as in Vandali 2011 for female and male speakers
    % Error rate, PB:
    [error_rate_PB parPB] = experiment_report_20140228_F0_as_Vandali(info_PB);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST-f Evaluation

info_LIST = info;
info_LIST = Ensure_field(info_LIST,'F0reference','praat');
info_LIST = Ensure_field(info_LIST,'F0max',400);
info_LIST = Ensure_field(info_LIST,'bCleanSpeech',1); % corrected plots on 9/10/2014


if bLISTf_SSN
    
    info_LIST.speaker = 'wdz'; % Female 
    %info.results_dir_name   = 'results-praat-CP810';
    info_LIST.figures_dir_name   = 'figures-praat-CP810';

    info = Ensure_field(info, 'figures_folder', '~/Documenten/LaTeX_Docs/paper/figures/');
    
    %info_LIST.results_dir_name   = 'output-SSN';
    %info_LIST.results_dir = [info_LIST.root_dir info.results_dir_name delim];
    
    % if info_LIST.bCleanSpeech == 1;
    %     info_LIST.praat       = [info_LIST.root_dir 'praat' delim];
    % else
    %     info_LIST.praat       = [info_LIST.root_dir 'praat-CP810' delim];
    % end
    
    info_LIST.bAssess = 0;
    info_LIST.bAnalyse = 1;
    info_LIST.F0max = 300;
    [info_LIST outLIST] = experiment_report_20140228_Physical_validation_LIST(info_LIST);

    % End of LIST-f evaluation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots added to paper (on 2014 04 02)

    [error_rate_LISTf parLIST] = experiment_report_20140228_F0_as_Vandali(info_LIST);

    snrf = [99 20 15 10 5 0 -5];
    error_rate_PB_1 = interp1(parPB.snr,error_rate_PB(1,:),snrf,'linear');
    error_rate_PB_2 = interp1(parPB.snr,error_rate_PB(2,:),snrf,'linear');
    error_rate_LISTf = interp1(parLIST.snr,error_rate_LISTf,snrf,'linear');
    error_rate_tot = [error_rate_PB_1; error_rate_PB_2; error_rate_LISTf];

    h(end+1) = Plot_errorrate_F0mod(error_rate_tot);
    filename = [info.figures_folder 'ac-error-rate'];

    if info_LIST.bSave
        Saveas(h(end), filename);
    end
    
end

h = [];
close all

n = min([size(tmpPBm.m1,1) size(outLIST.m1,1)]);

%% Figure 2
figure
stPlot.YLabel   = 'vErr [%]';
stPlot.yLim     = [0 100];
h(end+1) = Plot_measure([tmpPBm.m1(1:n) tmpPBf.m1(1:n) outLIST.m1(1:n)],[tmpPBm.s1(1:n) tmpPBf.s1(1:n) outLIST.s1(1:n)], '', stPlot);
stPlot = [];

figure
stPlot.YLabel   = 'uvErr [%]';
stPlot.yLim     = [0 30];
h(end+1) = Plot_measure([tmpPBm.m2(1:n) tmpPBf.m2(1:n) outLIST.m2(1:n)],[tmpPBm.s2(1:n) tmpPBf.s2(1:n) outLIST.s2(1:n)], '', stPlot);
stPlot = [];

figure
stPlot.YLabel   = 'gErr [%]';
stPlot.yLim     = [0 30];
h(end+1) = Plot_measure([tmpPBm.m3(1:n) tmpPBf.m3(1:n) outLIST.m3(1:n)],[tmpPBm.s3(1:n) tmpPBf.s3(1:n) outLIST.s3(1:n)], '', stPlot);
stPlot = [];

figure
stPlot.YLabel    = 'f0Dev [\Delta Hz]';
h(end+1) = Plot_measure([tmpPBm.m4(1:n) tmpPBf.m4(1:n) outLIST.m4(1:n)],[tmpPBm.s4(1:n) tmpPBf.s4(1:n) outLIST.s4(1:n)], '', stPlot);
stPlot = [];

speakers = 'Speech-materials';

%%
% Written values:
% snr = [99 20 10 5 0 -5];

nDecimals = 2;

disp('LIST-f, vErr, uvErr, gErr')
Round([outLIST.snr'   outLIST.m1*outLIST.t_voiced/outLIST.t_total ...
                outLIST.m2*outLIST.t_unvoiced/outLIST.t_total ...
                outLIST.m3*outLIST.t_voiced/outLIST.t_total],nDecimals)

% snrPB = [99 20 10 5 0 -5 -10 -15];
disp(' ')
disp('PB-f, vErr, uvErr, gErr')
Round([tmpPBf.snr'    tmpPBf.m1*tmpPBf.t_voiced/tmpPBf.t_total ...
                tmpPBf.m2*tmpPBf.t_unvoiced/tmpPBf.t_total ...
                tmpPBf.m3*tmpPBf.t_voiced/tmpPBf.t_total], nDecimals)

% snrPB = [99 20 10 5 0 -5 -10 -15];
disp(' ')
disp('PB-m, vErr, uvErr, gErr')
Round([tmpPBm.snr' tmpPBf.m1*tmpPBm.t_voiced/tmpPBm.t_total ...
                   tmpPBf.m2*tmpPBm.t_unvoiced/tmpPBm.t_total ...
                   tmpPBf.m3*tmpPBm.t_voiced/tmpPBm.t_total],nDecimals)

Avg1 = mean([tmpPBf.m1*tmpPBf.t_voiced/tmpPBf.t_total tmpPBm.m1*tmpPBm.t_voiced/tmpPBm.t_total]');
Avg2 = mean([tmpPBf.m2*tmpPBf.t_unvoiced/tmpPBf.t_total tmpPBm.m2*tmpPBm.t_unvoiced/tmpPBm.t_total]');
Avg3 = mean([tmpPBf.m3*tmpPBf.t_voiced/tmpPBf.t_total tmpPBm.m1*tmpPBm.t_voiced/tmpPBm.t_total]');

disp('Error rate PB female and male')
Round([tmpPBm.snr' error_rate_PB'],nDecimals)

disp(' ')
disp('PB-f+m, vErr, uvErr, gErr, total')
Round([tmpPBf.snr'  Avg1' ...
                    Avg2' ...
                    Avg3' ...
                    (Avg1+Avg2+Avg3)'], nDecimals)

if info.bSave
    for k = 1:length(h)
        try
            Saveas(h(k),[info.figures_folder speakers '-' num2str(k)]);
        catch
            Saveas(h(k),[speakers '-' num2str(k)]);
        end
    end
end

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
    stPlot.SeriesLabel  = {'PB male','PB female','LIST-f'};
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