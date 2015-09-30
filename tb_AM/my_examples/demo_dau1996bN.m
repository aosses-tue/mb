function outs = demo_dau1996bN(options)
% function outs = demo_dau1996bN(options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b. If stimuli are not
%       found, then they are along this script generated.
% 
% 2. Stand-alone example:
%       options.bSave = 1;
%       options.stim_durations = [10 20 40]; % ms
%       demo_dau1996b(options);
% 
%       options.bSave = 1;
%       options.stim_durations = [10 20 40 80 160 320]; % ms
%       demo_dau1996b(options);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 22/10/2014 % Update this date manually
% Last use on   : 13/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = Ensure_field(options, 'nExperiment',2); % Exp. 3 - signal integration
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1); % just main plot
options = Ensure_field(options, 'method','dau1996'); % dau1996  uses dau1996preproc
                                                     % dau1996a uses dau1996apreproc
nExperiment = options.nExperiment;
bPlot = 0; % This is for getting a lot of plots

options = Ensure_field(options, 'dB_SPL'      , 76);
options = Ensure_field(options, 'dB_SPL_noise', 77);
options = Ensure_field(options, 'output_dir', Get_TUe_paths('outputs'));
paths.outputs   = options.output_dir;

h = []; % we initialise handle for Figures

%% II.A Deterministic maskers: simultaneous masking

% Common parameters:
dB_SPL_noise    = options.dB_SPL_noise; % reference: Left audio file
dB_SPL_target   = options.dB_SPL_target; % target dB SPL
dB_SPL          = options.dB_SPL; 
dB_Gain         = dB_SPL_target - dB_SPL;

opts.method = options.method;

filenames = options.filenames;
switch nExperiment
          
	case 20
    %% 0. Backward masking in frozen-noise maskers
        
        % test stim - 10 ms, same stimuli than Forward masking (experiment IIB0)
        options     = Ensure_field(options, 'stim_durations',10);
        
        % [noise_onset, noise_dur, t_silence_aft, t_total_duration] = Create_noise_dau1996_default(nExperiment);
        
        % noise_offset        = noise_onset + noise_dur;
        test_dur            = options.stim_durations*1e-3;
        test_onsets         = options.test_onsets;
        N_conditions        = length(test_onsets);
        
        opts.fc_idx         = 1000;
        
        stim_durations      = options.stim_durations;
                
        [innoise fs] = Wavread(filenames{end});
        options.bSave_noise = 0;
        
        for j = 1:length(filenames)-1
            [instim(:,j) fs] = Wavread(filenames{j});
            instim(:,j) = gaindb(instim(:,j),dB_Gain);
        end
        options.fs = fs;
        
        opts.bPlot      = bPlot;
        
        for i = 1:N_conditions
            
            time_offset = max( 50e-3, abs(min(test_onsets(i))) ); % samples to be added to both noise and stim
            tmp_noise = Gen_silence(               time_offset,fs);
            N_added_noise   = length(tmp_noise);
            
            tmp_innoise = [tmp_noise; innoise(1:end-N_added_noise)];
            L = length(tmp_innoise);
            exp1 = sprintf('out_stim%.0f = Dau1996compare(tmp_innoise,instim(1:L,i),fs,opts);',i);
            exp2 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig1'');',i,i);
            exp3 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig2'');',i,i);
            exp4 = sprintf('outs.out_stim%.0f = out_stim%.0f;',i,i);
            eval(exp1);
            eval(exp2);
            eval(exp3);
            eval(exp4);
            
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
