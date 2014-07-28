function VoD_get_initial_time_measured
% function VoD_get_initial_time_measured
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on : 24/6/2014
% Last update: 27/6/2014 % Update this date manually
% Last used  : 27/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subdir_db_vod    = Get_TUe_subpaths('db_voice_of_dragon');

dir_calibrated_m = subdir_db_vod.dir_calibrated_m; 
dir_calibrated_p = subdir_db_vod.dir_calibrated_p; 
dir_meas_def     = subdir_db_vod.dir_meas_def;

misc             = Get_VoD_params(1); % input = 0, to not store 'period' results
T                = misc.Tmodel;
last_take        = misc.last_take;
modes2check      = 2:5;
bHPF             = 1;

if bHPF
    lblFilter = 'filt'; % Suffix for wav-files, in case they were filtered
else
    lblFilter = '';
end

close all

for mode = modes2check
    
    mode_idx = mode-1;
    modus   = num2str(mode_idx);
    take    = last_take(mode_idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measured wav files:
    field = '1';
    filename{1} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    field = '2';
    filename{2} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    filename{3} = [dir_meas_def     'modus ' modus '_v' num2str(take) '-3.wav'];
    
    [yfar  fs] = Wavread(filename{1});
    [ynear fs] = Wavread(filename{2});
    y4period   = Wavread(filename{3});
    
    t = ( 0:length(ynear )-1 )/fs;
    
    tstart = misc.Ti_coil(mode_idx);
    
    figure;
    plot(t,ynear,t,y4period); xlim([0 1]), hold on
    [yy ty stats] = rmsdb_sec(ynear,t,tstart,fs);
    legend('near','4period')
    title( sprintf('Mode %.0f, take %.0f, time found = %3f [s]',mode,take,stats.tmin) )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modelled wav files:
    
    field = '1';
    filename{4} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
    field = '2';
    filename{5} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];

    [yfarp  fs] = Wavread(filename{4});
    [ynearp fs] = Wavread(filename{5});
        
    t = ( 0:length(yfar )-1 )/fs;
    
    tstart = misc.Ti_coil(mode-1);
    
    figure;
    plot(t,ynearp); xlim([0 1]), hold on
    [yy ty stats] = rmsdb_sec(ynearp,t,tstart,fs);
    legend('near')
    title( sprintf('Predicted: Mode %.0f, take %.0f, time found = %3f [s]',mode,take,stats.tmin) )
    % do not consider the stats.min for these measured values. Corrections 
    % were introduced visually to fit local minima with peaks across the 
    % period
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end