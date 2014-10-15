function [errorrate_tot outs] = experiment_report_20140228_F0_as_Vandali(info)
% function [errorate_tot outs] = experiment_report_20140228_F0_as_Vandali(info)
%
% Reads results inside info.results_dir:
%   - [speakers{i} '-ErrRateSim.txt'  ]
% 
% Programmed by Alejandro Osses, ExpORL 2014
% Last update on: 09/10/2014
% Last use on   : 09/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin == 0
    localpaths = Get_paths([],1); % addpath('~/Documenten/MATLAB/MATLAB_svn/Utility/')
    info = [];
end

info = Ensure_field(info, 'bAssess' , 1); % To assess physical measures
info = Ensure_field(info, 'bAnalyse', 1);
info = Ensure_field(info, 'bPlot'   , 0);
 
% close all
info        = Ensure_field(info,'root_dir'      ,'~/Documenten/fda_eval/');
info        = Ensure_field(info,'F0reference'   ,'laryngeal');
info        = Ensure_field(info,'F0max'         , 400);
info        = Ensure_field(info, 'Fs'           , round(15659.375));
info        = Ensure_field(info, 'speaker'      , 'sb');
 
if strcmp(info.F0reference, 'laryngeal')
    info.bCleanSpeech   = 1;
    info.nSentences = 50; % 50 sentences in PB database
end
 
info        = Ensure_field(info,'nSentences',350); % 350 sentences in LIST-f database
N = info.nSentences;

p.CFG.Fs        = info.Fs;
p.F0mod.F0max   = info.F0max;
 
% info.praat      = [info.root_dir 'praat'   delim];

if strcmp(info.F0reference,'laryngeal')
    
    info = Ensure_field(info, 'results_dir_name', 'results');
    info = Ensure_field(info, 'figures_dir_name', 'figures');
    
end
 
try
    info    = Ensure_field(info, 'results_dir', [info.root_dir info.results_dir_name delim]);
catch
    info    = Ensure_field(info, 'results_dir', info.output); % Since october 2014
end

info    = Ensure_field(info, 'figures'    , [info.root_dir info.figures_dir_name delim]);

info    = Ensure_field(info,'bPlot',0);
info    = Ensure_field(info,'bSave', 1); % F0mod_validation_PB_database
   
if info.bAnalyse
    speakers = {'sb','rl'};
    gender   = {' (female)', ' (male)'};
    
    if strcmp(info.speaker,'LIST-f') | strcmp(info.speaker,'wdz') % both for LIST-f
        speakers = {info.speaker,'rl'};
        idx2analyse = 1;
    else
        idx2analyse = 1:2;
    end
    
    col_time    = 2;
    col_time_v  = 3;
    col_time_uv = 4;
    
    errorrate_tot = [];
    
    for i = idx2analyse
        h = [];
        % vErrSim = import_physical_measure([info.results_dir speakers{i} '-vErrSim.txt'  ], 2, N+1);
        errorrate       = import_physical_measure([info.results_dir speakers{i} '-ErrRateSim.txt'  ], 2, N+1, 9);
        database_info   = import_physical_measure([info.results_dir speakers{i} '-database-time-info.txt'  ], 2, N+1, 6);
   
        
        tot_time_v  = sum(database_info(:,col_time_v));
        tot_time_uv = sum(database_info(:,col_time_uv));
        tot_time    = sum(database_info(:,col_time  ));
       
        % Pooling Error rate:
        error_data = errorrate(:,2:end);
        vector_t_voiced = repmat( database_info(:,col_time_v),1,size(error_data,2));
        errortmp = sum(  errorrate(:,2:end).* vector_t_voiced )/(tot_time_v);
        errorrate_tot = [errorrate_tot; errortmp]; %  
        
    end
    
    outs.tot_time = tot_time;
    outs.tot_time_v = tot_time_v;
    outs.tot_time_uv = tot_time_uv;
    
    outs.snr = Get_snr(errorrate_tot);
    outs.errorrate_tot = errorrate_tot;
    
    if info.bPlot
        h = Plot_errorrate(errorrate_tot);
        if info.bSave
            filename = [info.figures 'ac-vs-eTone2011'];
            Saveas(h, filename);
        end
    end
    
end

disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('EOF')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function snr = Get_snr(errorrate)

if size(errorrate,2) == 8
    snr       = [99 20    10 5 0 -5 -10 -15];
elseif size(errorrate,2) == 6
    snr       = [99 20    10 5 0 -5];
end

function h = Plot_errorrate(errorrate,snr)

% The following files were taken from Vandali2011, Figure 5
%       - Panel (a): male
%       - Panel (b): female
snr_vandali         = [99    12  8    4    0    -4  -8  -12];
error_vandali_f     = [2.38   3  4    5    7    13  30   56]; % Std criterion
error_vandali_m     = [3.98   6  7    9   10.5  14  21   30]; % Std criterion
error_vandali_f_SRT = [1.64   2  2.2  3    4     8  19   43]; % SRT > 0.6
error_vandali_m_SRT = [1.36   2  2.2  2    2.3   3   5   12]; % SRT > 0.6

% ----

errorrate_f     = errorrate(1,:);
errorrate_m     = errorrate(2,:);
if nargin < 2
    snr         = [99 20 15 10 5 0 -5 -10 -15];
end
snr_i           = [99 20    10 5 0 -5 -10 -15];

error_vandali_f     = interp1(snr_vandali, error_vandali_f    , snr, 'linear');
error_vandali_m     = interp1(snr_vandali, error_vandali_m    , snr, 'linear');
error_vandali_f_SRT = interp1(snr_vandali, error_vandali_f_SRT, snr, 'linear');
error_vandali_m_SRT = interp1(snr_vandali, error_vandali_m_SRT, snr, 'linear');
errorrate_f     = interp1(snr_i, errorrate_f, snr, 'linear');
errorrate_m     = interp1(snr_i, errorrate_m, snr, 'linear');

st.YLim         = [0 100];
st.YTick        = [0:10:100];
st.XLim         = [1 8.2];
st.XTickLabel   = {'Q','20','15','10','5','0','-5','-10',' '};

clGray = [0.5 0.5 0.5];

figure
subplot(1,2,1)
plot(1:length(snr), errorrate_m        , 'ko-.'              ,'MarkerFaceColor','k'   ,'LineWidth',1.5), hold on, grid on
%plot(1:length(snr), error_vandali_m    , 's--','Color',clGray,'MarkerFaceColor',clGray,'LineWidth',1.5)
plot(1:length(snr), error_vandali_m    , 'k>--'              ,'MarkerFaceColor','w'   ,'LineWidth',1.5)
plot(1:length(snr), error_vandali_m_SRT, 'k>-'               ,'MarkerFaceColor','k'   ,'LineWidth',1.5)

hl = legend('F0mod', 'eTone, Std.', 'eTone SRT > 0.6');
grid on

title('Male speech material')
ylabel('Error Rate (%)')
xlabel('SNR (dB)')

h   = gcf;
ha  = gca;

set(ha,'XLim'       ,st.XLim        )
set(ha,'YLim'       ,st.YLim        )
set(ha,'XTickLabel' ,st.XTickLabel  )
set(ha,'YTick'      ,st.YTick       )
set(hl,'Location'   , 'NorthWest'   )

%%%%
subplot(1,2,2)
plot(1:length(snr), errorrate_f        , 'ko-.'              ,'MarkerFaceColor','k'   ,'LineWidth',1.5), hold on, grid on
% plot(1:length(snr), error_vandali_f    , 's--','Color',clGray,'MarkerFaceColor',clGray,'LineWidth',1.5)
plot(1:length(snr), error_vandali_f    , 'k>--'              ,'MarkerFaceColor','w'   ,'LineWidth',1.5)
plot(1:length(snr), error_vandali_f_SRT, 'k>-'               ,'MarkerFaceColor','k'   ,'LineWidth',1.5)
grid on

title('Female speech material')
xlabel('SNR (dB)')

h   = gcf;
ha  = gca;

set(ha,'XLim'       ,st.XLim        )
set(ha,'YLim'       ,st.YLim        )
set(ha,'XTickLabel' ,st.XTickLabel  )
set(ha,'YTick'      ,st.YTick       )

set(h, 'Position', [1 25 1366 637])
