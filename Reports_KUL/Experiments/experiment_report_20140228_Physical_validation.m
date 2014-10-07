function info = experiment_report_20140228_Physical_validation(info)
% function info = experiment_report_20140228_Physical_validation(info)
%
% Programmed by Alejandro Osses, ExpORL 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    localpaths = Get_paths([],1); % addpath('~/Documenten/MATLAB/MATLAB_svn/Utility/')
    info = [];
end

info = Ensure_field(info, 'bAssess', 1); % To assess physical measures
info = Ensure_field(info, 'bAnalyse', 1);

close all
info        = Ensure_field(info,'root_dir'      ,'~/Documenten/fda_eval/');
info        = Ensure_field(info,'F0reference'   ,'laryngeal');
info        = Ensure_field(info,'F0max'         , 400);
info        = Ensure_field(info, 'Fs'           , round(15659.375));
info        = Ensure_field(info, 'speaker'      , 'sb');

if ~isfield(info,'bCleanSpeech')
    
    if strcmp(info.F0reference, 'laryngeal')
        info.bCleanSpeech   = 1; 
    elseif strcmp(info.F0reference,'praat')
        % info.bClearSpeech   = input('Press 1 if you want to use clean speech for Praat-F0 estimation. Press 0 if you want to use CP810 simulated speech: ');
        info.bCleanSpeech   = 0; % if 0 then emulated CP810 wav files will be used
    end
    
end

p.CFG.Fs        = info.Fs;
p.F0mod.F0max   = info.F0max;

if strcmp(info.F0reference,'laryngeal')
    
    info = Ensure_field(info, 'results_dir_name', 'results');
    info = Ensure_field(info, 'figures_dir_name', 'figures');
    info.praat      = [info.root_dir 'praat'   delim]; % Not used
    
else
    if info.bCleanSpeech == 0
        
        info = Ensure_field(info, 'results_dir_name', 'results-praat-CP810');
        info = Ensure_field(info, 'figures_dir_name', 'figures-praat-CP810');
        info = Ensure_field(info, 'praat', [info.root_dir 'praat-CP810' delim]);
        
    else
        
        info = Ensure_field(info, 'results_dir_name', 'results-praat-clean');
        info = Ensure_field(info, 'figures_dir_name', 'figures-praat-clean');
        info = Ensure_field(info, 'praat', [info.root_dir 'praat'       delim]);
        
    end
end

info    = Ensure_field(info, 'results_dir', [info.root_dir info.results_dir_name delim]);
info    = Ensure_field(info, 'figures'    , [info.root_dir info.figures_dir_name delim]);

info    = Ensure_field(info,'bPlot',0);
info    = Ensure_field(info,'bSave', 0); % F0mod_validation_PB_database
info.tstep       = 5e-4;

if info.bAssess
    
    info
    disp(['  Directory: ' info.root_dir      ' must exist (with corresponding subdirectories)'])
    disp(['  Directory: ' info.results_dir   ' might be replaced'])
    disp(['  Directory: ' info.figures       ' might be replaced'])
    disp('  ')
    disp(Display(['Press any button to continue']));
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause()
    
    Mkdir(info.results_dir);
    Mkdir(info.figures);
    
    info = F0mod_validation_PB_database(info, p);
   
end

if strcmp(info.speaker,'sb')
    idx_to_do = 1;
elseif strcmp(info.speaker,'rl')
    idx_to_do = 2;
end

Colors = [0 0 0; 1 1 1];

if info.bAnalyse
    speakers = {'sb','rl'};
    gender   = {' (female)', ' (male)'};
    
    for i = idx_to_do
        h = [];
         vErrSim    = import_physical_measure([info.results_dir speakers{i} '-vErrSim.txt'  ], 2, 51);
        uvErrSim    = import_physical_measure([info.results_dir speakers{i} '-uvErrSim.txt' ], 2, 51);
         gErrSim    = import_physical_measure([info.results_dir speakers{i} '-gErrSim.txt'  ], 2, 51);
        f0DevErrSim = import_physical_measure([info.results_dir speakers{i} '-f0ErrSim.txt' ], 2, 51);
        globalErrSim = import_physical_measure([info.results_dir speakers{i} '-totErrSim.txt' ], 2, 51);

         vErrNMT    = import_physical_measure([info.results_dir speakers{i} '-vErrNMT.txt'  ], 2, 51);
        uvErrNMT    = import_physical_measure([info.results_dir speakers{i} '-uvErrNMT.txt' ], 2, 51);
         gErrNMT    = import_physical_measure([info.results_dir speakers{i} '-gErrNMT.txt'  ], 2, 51);
        f0DevErrNMT = import_physical_measure([info.results_dir speakers{i} '-f0ErrNMT.txt' ], 2, 51);
        globalErrNMT = import_physical_measure([info.results_dir speakers{i} '-totErrNMT.txt' ], 2, 51);

        % [m s] = prepare_barplot(vErrNMT(:,2:end), vErrSim(:,2:end));
        [m s] = prepare_barplot(vErrSim(:,2:end));
        figure
        stPlot.yLim         = [0 100];
        stPlot.YLabel       = 'vErr [%]';
        h(end+1) = Plot_measure(m,s,[speakers{i} gender{i}], stPlot, Colors(i,:));
        stPlot = [];
        
        % [m s] = prepare_barplot(uvErrNMT(:,2:end), uvErrSim(:,2:end));
        [m s] = prepare_barplot(uvErrSim(:,2:end));
        figure
        stPlot.YLabel       = 'uvErr [%]';
        h(end+1) = Plot_measure(m,s,[speakers{i} gender{i}], stPlot, Colors(i,:));
        stPlot = [];
        
        % [m s] = prepare_barplot( gErrNMT(:,2:end),  gErrSim(:,2:end));
        [m s] = prepare_barplot( gErrSim(:,2:end));
        figure
        stPlot.YLabel       = 'gErr [%]';
        h(end+1) = Plot_measure(m,s,[speakers{i} gender{i}], stPlot, Colors(i,:));
        stPlot = [];
        
        % [m s] = prepare_barplot(f0DevErrNMT(:,2:end), f0DevErrSim(:,2:end));
        [m s] = prepare_barplot(f0DevErrSim(:,2:end));
        figure
        stPlot.YLabel       = 'f0Dev [\Delta Hz]';
        h(end+1) = Plot_measure(m,s,[speakers{i} gender{i}], stPlot, Colors(i,:));
        stPlot = [];
        
%         errTotNMT = vErrNMT+uvErrNMT+gErrNMT+f0DevErrNMT;
%         errTotSim = vErrSim+uvErrSim+gErrSim+f0DevErrSim;
% 
%         [m s] = prepare_barplot(globalErrNMT(:,2:end), globalErrSim(:,2:end));
%         figure
%         stPlot.yLim         = [0 54];
%         stPlot.YLabel       = 'Global error [%]';
%         h(end+1) = Plot_measure(m,s,[speakers{i} ' ErrTot'], stPlot);
%         stPlot = [];
        
        if info.bSave
            for k = 1:length(h)
                Saveas(h(k),[info.figures speakers{i} '-' num2str(k)]);
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('EOF')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m s] = prepare_barplot(x, y, z)
% function [m s] = prepare_barplot(x, y, z)
%
% m - mean
% s - std
m = [];
s = [];

if exist('x','var')
    [m1    ,s1    ]     = Get_mean(x);
%     [m1_tot,s1_tot]     = Get_mean(m1');
    m = [m m1'];
    s = [s s1'];
end
if exist('y','var')
    [m2    ,s2    ]     = Get_mean(y);
%     [m2_tot,s2_tot]     = Get_mean(m2');
    m = [m m2'];
    s = [s s2'];
end

if exist('z','var')
    [m3    ,s3    ]     = Get_mean(z);
%     [m3_tot,s3_tot]     = Get_mean(m3');
    m = [m m3'];
    s = [s s3'];
end

function h = Plot_measure(m,s,title_str,stPlot, Colors)

if nargin < 4
    stPlot = [];
end

nTicks = 10;
stPlot.figPos       = [0 0 1024 300];
stPlot.xTickLabel   = {'Q','20','10','5','0','-5'};
stPlot.xTick        = 1:length(stPlot.xTickLabel);
stPlot.xLim         = [ 0   max(stPlot.xTick)+1];
stPlot              = Ensure_field(stPlot, 'yLim', [-2 20]);
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
    stPlot.SeriesLabel  = {'NMT','Sim'};
    if ~exist('Colors','var')
        Colors              = [1 1 1; 0 0 0];  
    end
end

if isfield(stPlot, 'SeriesLabel')
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
else
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors);
end

grid on, hold on
plot(stPlot.xLim,[50 50],'r--')

h = gcf;
set(h,'Position', stPlot.figPos);

Handle = gca;
set(Handle,'YTick',stPlot.yTick)
set(Handle,'XTick',stPlot.xTick)
set(Handle,'XLim',stPlot.xLim)
set(Handle,'YLim',stPlot.yLim)
set(Handle,'XTickLabel',stPlot.xTickLabel)

set(h, 'PaperPositionMode','auto')
