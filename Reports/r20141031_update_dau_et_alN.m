function r20141031_update_dau_et_alN(options)
% function r20141031_update_dau_et_alN(options)
%
% 1. Description:
%       Applies either the dau1996a or the dau1996 model. The dau1996a is the
%       model as published in dau1996a (see Mendeley) which does not include
%       overshoot limiting. The dau1996 includes overshoot limiting.
%       Overshoot limitation was one of the proposed improvements, as stated
%       in the paper, which is how the authors explained the differences between
%       measured and predicted masking thresholds in the backward  masking
%       coondition.
% 
% 2. Stand-alone example:
%       % 2.1 Example:
%       r20141031_update_dau_et_al;
%
%       % 2.2 Example:
%       opts.method = 'dau1996';
%       r20141031_update_dau_et_al(opts);
% 
%       % 2.3 Example:
%       opts.method = 'dau1996a';
%       r20141031_update_dau_et_al(opts);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 27/10/2014
% Last update on: 18/03/2015 
% Last use on   : 29/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Common parameters:
options.dB_SPL_noise = 77;
criterion_corr  = 6.5; % Arbitrary number

% options = Ensure_field(options,'method','dau1996a'); % no overshoot limit
options = Ensure_field(options,'method','dau1996'); % overshoot limit = 10
options = Ensure_field(options,'bExpIC0',1); 
options = Ensure_field(options,'calc_method',1); % 1 - as originally programmed
                                                 % 2 - other normalisation
calc_method = options.calc_method; % calc_method is used internally to use different normalisations

options.criterion_corr = criterion_corr;

bExpIC0 = options.bExpIC0; % Backward masking     L613 on 18/03/2015, tested OK

p = Get_date;
options.output_dir  = options.pathaudio;
output_dir          = options.pathaudio; % everything stored in the same directory

h = [];

%% Backward masking: 
if bExpIC0 == 1; 

Threshold = [];
mue = [];
label_experiment    = 'Dau1996b-ExpIC0'; disp(label_experiment);
label_figure        = 'Backward masking experiment';

options.nExperiment = 20;

noise_onset = 0e-3; 
noise_dur   = 200e-3;

SPL_test = [10 16:10:76];
test_onset_ref  = [-20:5:10]*1e-3;

options.test_onsets = noise_onset+test_onset_ref;

options.dB_SPL_target   = 85;
options.dB_SPL          = 76;
options.compute = 'template'; % 'all'

options.bAddNoise = 0;
options.sigma = 0; warning('temporal')
outs85  = il_demo_dau1996bN(options); 
idx     = outs85.out_stim1.idx;
fs      = outs85.fs;

test_onsets = options.test_onsets;

N_conditions = length(test_onsets);
for i = 1:N_conditions
    exp1 = sprintf('template%.0f = outs85.out_stim%.0f.template;',i,i);
    eval(exp1);
end

for i = 1:length(SPL_test)
    options.dB_SPL_target = SPL_test(i);
    options.compute = 'difference';

    options.bAddNoise = 1;
    options.sigma = 0.8;
    outstest = il_demo_dau1996bN(options); 

    % Decision
    
    for j = 1:N_conditions
        
        switch calc_method
            case 1
                exp1 = sprintf('curr_diff = outstest.out_stim%.0f.curr_diff;',j); 
                exp2 = sprintf('[mue(j,i) tmp] = optimaldetector(curr_diff,template%.0f);',j);
            case 2
                exp1 = sprintf('curr_diff = outstest.out_stim%.0f.curr_diff;',j); 
                exp2 = sprintf('[mue(j,i) tmp] = optimaldetector(curr_diff,template%.0f,fs);',j);
        end
        eval(exp1);
        eval(exp2);
        disp('')
        
        bDebug = 1;
        if bDebug
            figure; plot(curr_diff,'r'); hold on; plot(template1); plot(tmp,'k'); legend('current difference','template',sprintf('CC (scaled) with CCV = %.2f',1/fs * sum(tmp)))
        end
    end
    
    disp('')
end

%% Decision making:
ptmp = Get_date;
filename_mat = [output_dir 'mue-' ptmp.date2print '.mat'];
save(filename_mat,'mue','options');

Threshold = demo_dau1996b_decision(mue,SPL_test,criterion_corr);

% Saving figures
figure;
plot(test_onset_ref*1000,Threshold,'o--'), grid on
xlabel('Signal onset relative to masker onset [ms]')
ylabel('Masked threshold [dB]')
title(sprintf('%s. Criterion = %.1f',label_figure,criterion_corr))

%
filename_Sref = [output_dir label_experiment '_Thres'];
Saveas(gcf,filename_Sref);

Thres_labels = [test_onsets; Threshold];

filename_Sref = [output_dir label_experiment '_Thres'];
save(filename_Sref,'Threshold');
disp(['Variable saved as: ' filename_Sref '.mat']);

filename_Sref = [output_dir label_experiment '_Thres_labels'];
save(filename_Sref,'Thres_labels');
disp(['Variable saved as: ' filename_Sref '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_outs = il_demo_dau1996bN(options)
% function il_outs = il_demo_dau1996bN(options)

options = Ensure_field(options, 'nExperiment',2); % Exp. 3 - signal integration
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1); % just main plot
options = Ensure_field(options, 'method','dau1996'); % dau1996  uses dau1996preproc
                                                     % dau1996a uses dau1996apreproc
options = Ensure_field(options,'compute','all');

opts.compute = options.compute;
opts.calc_method = options.calc_method;
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
        
        for j = 1:length(filenames)-1
            [instim(:,j) fs] = Wavread(filenames{j});
            instim(:,j) = gaindb(instim(:,j),dB_Gain);
        end
        options.fs = fs;
        
        opts.bPlot      = bPlot;
        
        for j = 1:N_conditions
            
            time_offset = max( 50e-3, abs(min(test_onsets(j))) ); % samples to be added to both noise and stim
            tmp_noise = Gen_silence(               time_offset,fs);
            N_added_noise   = length(tmp_noise);
            
            tmp_innoise = [tmp_noise; innoise(1:end-N_added_noise)];
            L = length(tmp_innoise);
            opts.bAddNoise = options.bAddNoise;
            opts.sigma = options.sigma;
            tmp = Dau1996compare(tmp_innoise,instim(1:L,j),fs,opts,opts.compute);
            exp1 = sprintf('out_stim%.0f = tmp;',j);
            eval(exp1);
            
            exp2 = sprintf('il_outs.out_stim%.0f = out_stim%.0f;',j,j);
            eval(exp2);
            
        end
        
        il_outs.fs = fs;
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%