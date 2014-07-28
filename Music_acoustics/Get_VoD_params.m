function params = Get_VoD_params(bGetT,bSave)
% function params = Get_VoD_params(bGetT,bSave)
%
% 1. Description:
%
%       params.last_take
%       params.last_take_idx
%       params.mf
%       params.t - time in seconds considered for the measured audio files
%       params.ti - NaN means that initial time has not been determined yet
%       params.ti_measured
%       params.ti_model
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       params = Get_VoD_params;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 18/06/2014
% Last update on: 28/07/2014 % Update this date manually
% Last used on  : 28/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bGetT = input('Do you want to determine the measured signal periods? (1 = yes; 0 = no): ');
end

if nargin < 2
    bSave = 0;
end

hFig = [];

if bSave 
    disp([mfilename '.m: you are going to save some results, if you want to cancel this press ctrl+C...'])
    disp('Continue?...')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause()
    disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.R      = 67/100; % rotation radio [m]
params.Xo_nf  = 83/100+params.R;
params.Yo_nf  = 0;
params.Zo_nf  = 204/100;

params.Xo_ff  = 256/100+params.R;
params.Yo_ff  = 0;
params.Zo_ff  = 168/100;

% Pos_obs = [Xo, Yo, Zo]; % Observer's position (microphone)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resonant frequencies as used in Hirschberg2013 (manually copied):
params.mf   = [     424.4, ... 
                    638.7, ...
                    851.8, ...
                    1061.7];

% Rotation periods as used in Hirschberg2013 (manually copied):
params.Tmodel = [   0.6023; ... % acoustic mode n = 2
                    0.3934; ... % acoustic mode n = 3
                    0.2961; ... % acoustic mode n = 4 % 0.2844
                    0.2688 ];   % acoustic mode n = 5
% NOTE: little difference for mode n = 2 for measured and modelled values. 
%       'Predicted' value (theory) was used as reference

params.last_take     = [2 3 4 3]; % number of takes for modus 1,2,3,4
params.last_take_idx = cumsum(params.last_take); % index if we pool the measurements
params.fs            = 10000;
params.t             = 0: 1/params.fs :20 - 2/params.fs; % 1 sample deleted

% NOTE: 2,3,4 and 3 takes respectively were available for acoustic modes 2 to 5.
%       I assumed that he last take was always the 'best'. This can be confirmed 
%       with the following analysis where you can see that measured signals 
%       corresponding to these takes are closer to theoretical values...

tmp = VoD_get_period_init;
params.ti_model     = tmp.ti_model;
params.ti_measured  = tmp.ti_measured;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated on: 01/07/2014
% Measured periods closer to global average (aligned near-field signals):
% (See "D:\MATLAB\Output\logs-MATLAB\log-VoD_comparisons-2014-7-1-at-15h-9m-48s")
params.meas_n2 = [10 11 12]; % cycle 12 added just to have 3 cycles
params.meas_n3 = [46 47 48 49]; % cycle 47, greater than 1, put here for continuity
params.meas_n4 = [61 62 63 64];
params.meas_n5 = [18 19 20 21 22];

params.nTmeas_stable = [ length(params.meas_n2); ... 
                         length(params.meas_n3); ...
                         length(params.meas_n4); ...
                         length(params.meas_n5) ];
params.tmeas_stable = params.nTmeas_stable .* params.Tmodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bGetT == 1
    
    subdir_db_vod = Get_TUe_subpaths('db_voice_of_dragon');
    
    sourcedir   = subdir_db_vod.dir_meas_def;
    filename    = Get_filenames(sourcedir,'*3.wav'); % Signals used to get periods
    T           = nan(length(filename),100); % 100 columns arbitrarily
    t_coil_all      = nan(length(filename),100); % 100 columns arbitrarily
    idx_count   = 4;
    for i = length(filename):-1:1
        
        [Counter misc] = Count_times_above_thr([sourcedir filename{i}]);
        
        T(i,1:length(misc.T)) = misc.T;
        t_coil_all(i,1:length(misc.ti)) = misc.ti;
        
        if i == params.last_take_idx(idx_count)
            Ti_coil(idx_count) = misc.ti(1);
            t_coil(idx_count,:) = t_coil_all(i,:);
            
            idx_count = idx_count - 1;
            if idx_count == 0
                idx_count = 1;
            end
        end
        
    end
    
    Perc = 70.8559/100;
    t_coil(:,1) = []; % removing incomplete cycles
    tstart = t_coil(:,1) - Perc*params.Tmodel;
    % params.ti_measured = tstart'; % uncomment to assign. These values were 
                                    % manually copied into VoD_get_period_init
    
    T = Delete_NaN_columns(T);
    t_coil_all = Delete_NaN_columns(t_coil_all);
    t_coil = Delete_NaN_columns(t_coil);
    
    N = Count_isnan(T);
    [m,s,ci] = Get_mean(T');
    
    params.T            = T;
    params.Ti_coil      = Ti_coil;
    params.t_coil_all   = t_coil_all;
    params.t_coil       = t_coil;
    
    MeasNr = (1:length(N))';
    
    % Information to be displayed
    params.info = [MeasNr N m' minmax(T)];
    
    % Real boxplot:
    [h   stats] = Boxplot(T');
    ylabel('Period [s]')
    xlabel('Measurement Number')
    
    hFig(end+1) = gcf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Uncomment the following lines to carry out the same analysis as presented
    % % In report sent on 20/06/2014
    %
    % TT = T - repmat(stats.Median',1,size(T,2)); 
    % figure;
    % Boxplot(TT');
    % ylabel('Period - Median(Period) [s]')
    % xlabel('Measurement Number')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Ration as discussed on 23/07/2014
    TT = T ./ repmat(stats.Median',1,size(T,2)); 
    figure;
    Boxplot(TT');
    ylabel('Period / Median(Period) [dimensionless]')
    xlabel('Measurement Number')
    
    hFig(end+1) = gcf;
    
    if bSave
        disp('Copy the following into your LaTeX file: ')
        latex(params.info);
        
        outputpath = Get_TUe_paths('outputs');
        for i = 1:length(hFig)
            Saveas(hFig(i),[outputpath mfilename '-' num2str(i)])
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end