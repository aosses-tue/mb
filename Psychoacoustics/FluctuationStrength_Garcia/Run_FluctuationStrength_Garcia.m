function Run_FluctuationStrength_Garcia
% function Run_FluctuationStrength_Garcia
%
% 1. Description:
%
% 2. Stand-alone example:
%       Run_FluctuationStrength_Garcia;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 03/02/2016
% Last update on: 03/02/2016 
% Last use on   : 03/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
N  = 1*44100; % 1 sec at 44100 Hz

Run_mode = 0; % 0 = as ran by Rodrigo

bAM_fmod = 1;
bFM_fmod = 0;
dirout = 'd:\Documenten-TUe\02-Experiments\2015-APEX-Rodrigo\APEX_shared\Stimuli\';

if bAM_fmod
    
    tests = [0 0.5 1 2 4 8 16];
    fnames = {  [dirout 'AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']};
%     fnames = {  [dirout 'AM_tone-fm_0.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']
%                 [dirout 'AM_tone-fm_0.50-fc_1000-md_40-SPL_70-w_25-fs_44100-N_176400.wav']; ...
%                 [dirout 'AM_tone-fm_1.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_132300.wav']; ...
%                 [dirout 'AM_tone-fm_2.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ...
%                 [dirout 'AM_tone-fm_4.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ...
%                 [dirout 'AM_tone-fm_8.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']; ...
%                 [dirout 'AM_tone-fm_16.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_88200.wav']};

    switch Run_mode
        case 0
            % delete([Get_path('FluctuationStrength_Model_Data') '*.mat'])
            params = Get_fluctuation_strength_params(N,fs);
            params.debug = 'none';
        case 1
            load('params-20160203-1914.mat');
            load('specs-20160203-1917.mat'); % specs = 1x8 cell
    end

    % % Imported from function Get_specs:
    % % handles = {@spec_AM_fmod @AM_fc @AM_md @AM_SPL @FM_fm @FM_fc @FM_df @FM_SPL};
    % handles = {@spec_AM_fmod};
    % specs   = cellfun(@(x)x(),handles,'UniformOutput',false);

    tic;
    model_AM_fmod = cellfun(@(s)il_process_spec(params,s),fnames);
    toc;
    
    figure;
    plot( 1:length(tests),model_AM_fmod ); grid on
    set(gca,'XTick',1:length(tests))
    set(gca,'XTickLabel',tests);
end

if bFM_fmod
    
    tests = [0 0.5 1 2 4 8 16];
    fnames = {  [dirout 'FM_tone-fm_0.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']
                [dirout 'FM_tone-fm_0.50-fc_1500-df_700-SPL_70-w_25-fs_44100-N_176400.wav']; ...
                [dirout 'FM_tone-fm_1.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_132300.wav']; ...
                [dirout 'FM_tone-fm_2.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
                [dirout 'FM_tone-fm_4.00-fc_1500-df_700-SPL_40-w_25-fs_44100-N_88200.wav']; ...
                [dirout 'FM_tone-fm_8.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']; ...
                [dirout 'FM_tone-fm_16.00-fc_1500-df_700-SPL_70-w_25-fs_44100-N_88200.wav']};

    switch Run_mode
        case 0
            % delete([Get_path('FluctuationStrength_Model_Data') '*.mat'])
            params = Get_fluctuation_strength_params(N,fs);
            params.debug = 'none';
        case 1
            load('params-20160203-1914.mat');
            load('specs-20160203-1917.mat'); % specs = 1x8 cell
    end

    % % Imported from function Get_specs:
    % % handles = {@spec_AM_fmod @AM_fc @AM_md @AM_SPL @FM_fm @FM_fc @FM_df @FM_SPL};
    % handles = {@spec_AM_fmod};
    % specs   = cellfun(@(x)x(),handles,'UniformOutput',false);

    tic;
    model_AM_fmod = cellfun(@(s)il_process_spec(params,s),fnames);
    toc;
    
    figure;
    plot( 1:length(tests),model_AM_fmod ); grid on
    set(gca,'XTick',1:length(tests))
    set(gca,'XTickLabel',tests);
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = il_process_spec(params,fnames)
    
    % % From Get_model_data:
    % fnames = Get_stimuli_filenames(specs,params.N);
    model = Get_fluctuation_strength(fnames,params);
    
    disp('')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
