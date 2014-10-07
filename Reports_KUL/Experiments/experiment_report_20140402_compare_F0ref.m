function experiment_report_20140402_compare_F0ref(info)
% function experiment_report_20140402_compare_F0ref(info)
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPB_male = 0;
bPB_female = 0;
bLISTf_white = 0;
bLISTf_SSN = 0;

info.bAssess    = 0;
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

info.F0reference        = 'laryngeal';
info.F0max              = 400;

info.bCleanSpeech       = 1; % if 0: CP810 simulations will be used
info.bPlot              = 0;
info.speaker            = 'sb'; % Female 
info.results_dir_name   = 'results';
info.figures_dir_name   = 'figures';


Colors = [0 0 0; 1 1 1];
h = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PB database
if bPB_female
    [info tmp2] = experiment_report_20140402_Physical_validation(info);
    close all
end

if bPB_male
    info.speaker            = 'rl'; % Male

    [info tmp1] = experiment_report_20140402_Physical_validation(info);
    close all
end

if bPB_male & bPB_female
    % Analysis as in Vandali 2011 for female and male speakers
    % Error rate, PB:
    error_rate_PB = experiment_report_20140228_F0_as_Vandali(info);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST-f Evaluation

info.F0reference        = 'praat';
info.F0max              = 400;
info.bCleanSpeech       = 0;

if bLISTf_SSN
    
    info.speaker            = 'sb'; % Female 
    info.results_dir_name   = 'results-praat-CP810';
    info.figures_dir_name   = 'figures-praat-CP810';

    info = Ensure_field(info, 'figures_folder', '~/Documenten/LaTeX_Docs/paper/figures/');

    info_LIST = info;
    info_LIST.bAssess = 0;
    info_LIST.bAnalyse = 1;
    info_LIST.F0max = 300;
    [info_LIST tmp3] = experiment_report_20140228_Physical_validation_LIST(info_LIST);

    % End of LIST-f evaluation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots added to paper (on 2014 04 02)

    error_rate_LISTf = experiment_report_20140228_F0_as_Vandali(info_LIST);

    error_rate_tot = [error_rate_PB; error_rate_LISTf];

    h(end+1) = Plot_errorrate_F0mod(error_rate_tot);
    filename = [info.figures_folder 'ac-error-rate'];

    if info_LIST.bSave
        Saveas(h(end), filename);
    end
    
end

h = [];
close all

figure
stPlot.YLabel   = 'vErr [%]';
stPlot.yLim     = [0 100];
h(end+1) = Plot_measure([tmp1.m1 tmp2.m1 tmp3.m1],[tmp1.s1 tmp2.s1 tmp3.s1], '', stPlot);
stPlot = [];

figure
stPlot.YLabel   = 'uvErr [%]';
h(end+1) = Plot_measure([tmp1.m2 tmp2.m2 tmp3.m2],[tmp1.s2 tmp2.s2 tmp3.s2], '', stPlot);
stPlot = [];

figure
stPlot.YLabel    = 'gErr [%]';
h(end+1) = Plot_measure([tmp1.m3 tmp2.m3 tmp3.m3],[tmp1.s3 tmp2.s3 tmp3.s3], '', stPlot);
stPlot = [];

figure
stPlot.YLabel    = 'f0Dev [\Delta Hz]';
h(end+1) = Plot_measure([tmp1.m4 tmp2.m4 tmp3.m4],[tmp1.s4 tmp2.s4 tmp3.s4], '', stPlot);
stPlot = [];

speakers = 'Speech-materials';

if info.bSave
    for k = 1:length(h)
        Saveas(h(k),[info.figures speakers '-' num2str(k)]);
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