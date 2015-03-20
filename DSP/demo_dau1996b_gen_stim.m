function outs = demo_dau1996b_gen_stim(nExperiment,options)
% function outs = demo_dau1996b_gen_stim(nExperiment,options)
%
% 1. Description:
%
% 2. Stand-alone example:
%       demo_dau1996b_gen_stim;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/03/2015
% Last update on: 13/03/2015 % Update this date manually
% Last use on   : 18/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    options = [];
end

if nargin < 1
    nExperiment = 2;
end

options = Ensure_field(options,'bSave_noise' ,1);
options = Ensure_field(options,'bSave'       ,1);
options = Ensure_field(options,'dB_SPL_noise',77);
options = Ensure_field(options,'dB_SPL'      ,85);
options = Ensure_field(options,'fs'          ,44100);
options = Ensure_field(options,'output_dir'  ,Get_TUe_paths('outputs'));

N_stim = 6;

bSave_noise = options.bSave_noise;
bSave       = options.bSave;
output_dir  = options.output_dir;

bDiary = bSave;
Diary(mfilename,bDiary);

dB_SPL      = options.dB_SPL;
fs          = options.fs;

[noise_onset, t_duration, t_silence_aft, t_total_duration] = Create_noise_dau1996_default(nExperiment);

%% Noise generation/or wav read:
switch nExperiment % Noise generation
    case {1 2 3}
        filename    = [output_dir 'dau1996b_expI_noisemasker']; % used in Exp 2, 3
        
    case 20
        filename    = [output_dir 'dau1996b_expIB0_noisemasker']; % used in Exp 20    
end

if bSave_noise == 1
    innoise = Create_noise_dau1996(nExperiment,filename,options);
else
    innoise = Wavread(filename);
end

outs.innoise = innoise;

%%
switch nExperiment
    
    case 2 % Frozen noise
        
        % Generating the test tones:
        %   Common stim params
        onset   = noise_onset + 115e-3; % temporarily onset equal to the one of the noise + 115 ms
        f       = 1000;
        win     = 1; % 1 = Hanning window
        dur     = 5e-3;
        options.test_phases = [0:2/8:2]; % times pi
        test_phases = options.test_phases;
        N_conditions = length(test_phases);
        
        filename1   = [output_dir 'dau1996b_expI2_stim-5ms-' num2str(options.dB_SPL) '-phase-0-pi'];
        for i = 2:N_conditions
            exp1 = sprintf('filename%.0f = [output_dir ''dau1996b_expI2_stim-5ms-'' num2str(options.dB_SPL) ''-phase-%.0f_10-pi''];',i,test_phases(i)*10);
            eval(exp1);
        end
        
        % Stim 1
        [instim1, t1] = Create_sin4this_exp(f,test_phases(1)*pi,dur,fs,win,onset,dB_SPL,t_total_duration);
        if bSave == 1
            Wavwrite(instim1,fs,filename1); % Save instim1
        end

        for i = 2:N_conditions
            exp1 = sprintf('[instim%.0f, t%.0f] = Create_sin4this_exp(f,%.4f*pi,dur,fs,win,onset,dB_SPL,t_total_duration);',i,i,test_phases(i));
            exp2 = sprintf('Wavwrite(instim%.0f,fs,filename%.0f);',i,i);
            eval(exp1);
            if bSave == 1
                eval(exp2); % Save instim2-instim6
            end
        end
        
        outs.instim1 = instim1;
        outs.instim2 = instim2;
        outs.instim3 = instim3;
        outs.instim4 = instim4;
        outs.instim5 = instim5;
        outs.instim6 = instim6;
        outs.instim7 = instim7;
        outs.instim8 = instim8;
        outs.instim9 = instim9;
        
    case 3 % Frozen noise, signal fc=3 kHz, onset at 100 ms
        
        options.stim_durations = [10 20 40 80 160 320]; % ms
        stim_durations = options.stim_durations;
        
        filename1         = [output_dir 'dau1996b_expI3_stim01-' num2str(options.dB_SPL)];
        filename2         = [output_dir 'dau1996b_expI3_stim02-' num2str(options.dB_SPL)];
        filename3         = [output_dir 'dau1996b_expI3_stim03-' num2str(options.dB_SPL)];

        if N_stim > 3
            filename4     = [output_dir 'dau1996b_expI3_stim04-' num2str(options.dB_SPL)];
            filename5     = [output_dir 'dau1996b_expI3_stim05-' num2str(options.dB_SPL)];
            filename6     = [output_dir 'dau1996b_expI3_stim06-' num2str(options.dB_SPL)];
        end
        
        % Generating the test tones:
        %   Common stim params
        f       = 3000;
        onset   = noise_onset + 100e-3; % noise_onset is 0 by default
        win     = 1; % 1 = Hanning window

        % Stim 1
        dur     = 10e-3;
        [instim1, t1] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
        if bSave == 1
            Wavwrite(instim1,fs,filename1);
        end

        % Stim 2
        dur     = 20e-3;
        [instim2, t2] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
        if bSave == 1
            Wavwrite(instim2,fs,filename2);
        end

        % Stim 3
        dur     = 40e-3;
        [instim3, t3] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
        if bSave == 1
            Wavwrite(instim3,fs,filename3);
        end

        if N_stim > 3
            
            dur     = stim_durations(4)*1e-3;
            [instim4, t4] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            if bSave == 1
                Wavwrite(instim4,fs,filename4);
            end

            dur     = stim_durations(5)*1e-3;
            [instim5, t5] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            if bSave == 1
                Wavwrite(instim5,fs,filename5);
            end

            dur     = stim_durations(6)*1e-3;
            [instim6, t6] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            if bSave == 1
                Wavwrite(instim6,fs,filename6);
            end
            
            outs.instim4 = instim4;
            outs.instim5 = instim5;
            outs.instim6 = instim6;
            
        end
        
        outs.instim1 = instim1;
        outs.instim2 = instim2;
        outs.instim3 = instim3;
        
    case 20
        
        %% Generating the test tones
        f           = 1000;
        win         = 1; % 1 = Hanning window
        test_dur    = 10e-3;
        test_onsets = noise_onset+ [-100:40:-20,-15:5:15]*1e-3; %[-100:30:-40,-20:4:20]*1e-3;
        N_conditions = length(test_onsets);
        
        onset       = noise_onset; % temporarily onset equal to the one of the noise
        
        filename0   = [output_dir 'dau1996b_expIB0_stim-10ms-' num2str(options.dB_SPL)];
        filename1   = [output_dir 'dau1996b_expIB0_stim-10ms-' num2str(options.dB_SPL) '-1'];
        
        for i = 2:N_conditions
            exp0 = sprintf('filename%.0f   = [output_dir ''dau1996b_expIB0_stim-10ms-'' num2str(options.dB_SPL) ''-%.0f''];',i,i);
            eval(exp0)
        end
        
        [instim0, t0] = Create_sin4this_exp(f,0,test_dur,fs,win,onset,dB_SPL,t_total_duration);
        
        if length(instim0) ~= length(innoise);
            instim0 = [instim0; zeros(length(innoise)-length(instim0),1)];
        end
        
        if bSave == 1
            Wavwrite(instim0,fs,filename0);
        end
        
        for i = 1:N_conditions
            
            time_offset = max( 50e-3, abs(min(test_onsets(i))) ); % samples to be added to both noise and stim
            tmp_noise = Gen_silence(               time_offset,fs);
            tmp_insig = Gen_silence(test_onsets(i)+time_offset,fs);
            N_added         = length(tmp_insig); 
            N_added_noise   = length(tmp_noise);
            
            tmp_innoise = [tmp_noise; innoise(1:end-N_added_noise)];
            exp0 = sprintf('instim%.0f = [tmp_insig; instim0(1:end-N_added)];',i);
            exp1 = sprintf('instim%.0f = [tmp_insig; instim0(1:end-N_added)];',i);
            eval(exp0);
            eval(exp1);
            
            if bSave == 1
                exp2 = sprintf('Wavwrite(instim%.0f,fs,filename%.0f);',i,i);
                eval(exp2)
            end
            exp3 = sprintf('outs.instim%.0f = instim%.0f;',i,i);
            eval(exp3);            
            
        end
            
end
            
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [y,t] = Create_sin4this_exp(f,start_phase,dur,fs,win,onset,SPL,total_duration)

    [y, t]= Create_sin_phase(f,start_phase,dur,fs,win);
    
    % (f,dur,fs,win,onset,dB_SPL_above_thr);
    y  = setdbspl(y,SPL);
    y  = [Gen_silence(onset,fs); y]; 
    
    try % Append silence only if total_duration has been specified
    	y = [y; Gen_silence(total_duration-max(t)-onset-1/fs,fs)];
    end

    t = (1:length(y))/fs; % redefine t
    
