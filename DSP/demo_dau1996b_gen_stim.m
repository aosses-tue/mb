function outs = demo_dau1996b_gen_stim(nExperiment,options)
% function outs = demo_dau1996b_gen_stim(nExperiment,options)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/03/2015
% Last update on: 13/03/2015 % Update this date manually
% Last use on   : 13/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    options = [];
end

if nargin < 1
    nExperiment = 3;
end

options = Ensure_field(options,'bSave_noise' ,1);
options = Ensure_field(options,'dB_SPL_noise',77);
options = Ensure_field(options,'dB_SPL'      ,85);
options = Ensure_field(options,'fs'          ,44100);

N_stim = 6;

bSave       = options.bSave_noise;

bDiary = bSave;
Diary(mfilename,bDiary);

dB_SPL      = options.dB_SPL;
fs          = options.fs;

switch nExperiment
    case 3
        
        options.stim_durations = [10 20 40 80 160 320]; % ms
        stim_durations = options.stim_durations;
        
        filename          = 'dau1996b_expI3_noisemasker'; % used in Exp 3
        filename1         = ['dau1996b_expI3_stim01-' num2str(options.dB_SPL)];
        filename2         = ['dau1996b_expI3_stim02-' num2str(options.dB_SPL)];
        filename3         = ['dau1996b_expI3_stim03-' num2str(options.dB_SPL)];

        if N_stim > 3
            filename4     = ['dau1996b_expI3_stim04-' num2str(options.dB_SPL)];
            filename5     = ['dau1996b_expI3_stim05-' num2str(options.dB_SPL)];
            filename6     = ['dau1996b_expI3_stim06-' num2str(options.dB_SPL)];
        end
        
        [noise_onset, t_duration, t_silence_aft, t_total_duration] = Create_noise_dau1996_default(nExperiment);

        if bSave == 1
            innoise = Create_noise_dau1996(nExperiment,filename,options);
        end

        % Generating the test tones
        % Common stim params
        onset   = noise_onset + 100e-3; % noise_onset is 0 by default
        f       = 3000;
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
        
        outs.innoise = innoise;
        outs.instim1 = instim1;
        outs.instim2 = instim2;
        outs.instim3 = instim3;
        
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
    
