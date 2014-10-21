function outs = demo_dau1996b(options)
% function outs = demo_dau1996b(options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 1;
%       options.stim_durations = [10 20 40]; % ms
%       demo_dau1996b(options);
% 
%       options.bSave = 1;
%       options.stim_durations = [10 20 40 80 160 320]; % ms
%       demo_dau1996b(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 20/10/2014 % Update this date manually
% Last use on   : 20/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO: 4. Signal frequency

if nargin == 0
    close all
    options = [];
end

options = Ensure_field(options, 'nExperiment',3); % Exp. 3 - signal integration
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1); % just main plot

nExperiment = options.nExperiment;
bPlot = 0; % This is for getting a lot of plots

options = Ensure_field(options, 'dB_SPL'      , 85);
options = Ensure_field(options, 'dB_SPL_noise', 77);

h = []; % we initialise handle for Figures
paths.outputs   = Get_TUe_paths('outputs');
        
%% II.A Deterministic maskers: simultaneous masking

% Common parameters:
dB_SPL_noise    = options.dB_SPL_noise; % reference: Left audio file
dB_SPL          = options.dB_SPL;

switch nExperiment
       
    case 1
    %% 1. Temporal position of short signals in frozen-noise maskers
        options     = Ensure_field(options, 'stim_durations',5);
        options     = Ensure_field(options,'test_onsets',[95:10:195]*1e-3);
        test_onsets         = options.test_onsets;
        opts.fc_idx         = 1000;
        
        stim_durations      = options.stim_durations;
        N_stim              = length(stim_durations);
                
        infilename          = 'dau1996b_expII1_noisemasker'; % used in Exp 1 and 2, 10
        infilename0        = ['dau1996b_expII1_stim-5ms-' num2str(options.dB_SPL)];
        
        try 
            [innoise fs] = Wavread([paths.outputs infilename  '.wav']);
        catch
            options.bGenerate   = 1;
            options.bSave_noise = 1;
        end
        
        try
            [instim0 fs] = Wavread([paths.outputs infilename0 '.wav']);
            options.bGenerate   = 0;
        catch
            options.bGenerate   = 1;
            fs          = 48000;
            options     = Ensure_field(options,'bSave_noise',0);
        end
        
        options.fs = fs;
        
        filename    = [paths.outputs infilename];
        filename0   = [paths.outputs infilename0];
        
        if options.bGenerate
        
            [t_silence_bef, t_duration, t_silence_aft, t_total_duration] = Create_noise_default(tagNoise);
            
            if options.bSave_noise == 1
                tagNoise = 1;
                Create_noise(tagNoise,filename,options);
            end
            
            %% Generating the test tones

            onset   = t_silence_bef; % temporarily onset equal to the one of the noise
            f       = 1000;
            win     = 1; % 1 = Hanning window

            % Stim 1
            dur     = 5e-3;
            [instim0, t0] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            Wavwrite(instim0,fs,filename0);

        end

        opts.bPlot      = bPlot;
        
        N_conditions = length(test_onsets);
        
        for i = 1:N_conditions
            
            tmp_insig = Gen_silence(test_onsets(i),fs);
            N_added = length(tmp_insig); 
            
            exp1 = sprintf('instim%.0f = [tmp_insig; instim0(1:end-N_added)];',i);
            exp2 = sprintf('out_stim%.0f = Dau1996compare(innoise,instim%.0f,fs,opts);',i,i);
            
            % idx = out_stim1.idx;
            
            eval(exp1);
            eval(exp2);
            
            if N_conditions > 3 % remove fields to avoid 'out of memory'
                
                exp3 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig1'');',i,i);
                exp4 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig2'');',i,i);
                eval(exp3);
                eval(exp4);
                
            end
            
            exp5 = sprintf('outs.out_stim%.0f = out_stim%.0f;',i,i);
            eval(exp5);
            
        end
        
    case 2
    %% 2. Relative phase
        % 1kHz tone, onset = 115 ms, phase_i between 0 and 2*pi
        infilename          = 'dau1996b_expII1_noisemasker'; % used in Exp 1 and 2
        infilename1        = ['dau1996b_expII2_stim-5ms-' num2str(options.dB_SPL) '-phase-0-pi'];
        opts.fc_idx         = 1000;
        
        options = Ensure_field(options,'test_phases',[0:2/8:2]);
        test_phases = options.test_phases;
        
        N_conditions = length(test_phases);
        for i = 2:N_conditions
            exp1 = sprintf('infilename%.0f = [''dau1996b_expII2_stim-5ms-'' num2str(options.dB_SPL) ''-phase-%.0f_10-pi''];',i,test_phases(i)*10);
            eval(exp1);
        end
           
        try 
            [innoise fs] = Wavread([paths.outputs infilename  '.wav']);
        catch
            options.bGenerate   = 1;
            options.bSave_noise = 1;
        end
        
        try
            [instim1 fs] = Wavread([paths.outputs infilename1 '.wav']);
            for i = 2:N_conditions
                exp1 = sprintf('[instim%.0f,fs] = Wavread([paths.outputs infilename%0.f ''.wav'']);',i,i);
                eval(exp1);
            end
            options.bGenerate   = 0;
        catch
            options.bGenerate   = 1;
            options     = Ensure_field(options,'bSave_noise',0);
            fs          = 48000;
        end
        
        options.fs = fs;
        
        filename    = [paths.outputs infilename];
        filename1   = [paths.outputs infilename1];
        for i = 2:N_conditions
            exp1 = sprintf('filename%.0f = [paths.outputs infilename%0.f];',i,i);
            eval(exp1);
        end
        if options.bGenerate
            
            [t_silence_bef, t_duration, t_silence_aft, t_total_duration] = Create_noise_default(tagNoise);
            
            if options.bSave_noise == 1
                tagNoise = 2;
                Create_noise(tagNoise,filename,options);
            end
            
            % Nsil_bef    = round(options.fs*t_silence_bef);
            % Nnoise      = round(options.fs*t_duration);
            % Nsil_aft    = round(options.fs*t_silence_aft);

            % if options.bSave_noise == 1
            %     % Gen1: white noise, band-pass filtered
            %     y = wgn(Nnoise,1,1);
            % 
            %     y   =  y(:); % ensures it is a column vector
            % 
            %     Wn = [20 5000]/(options.fs/2); % Normalised cutoff frequency        
            %     [b,a] = butter(4,Wn); % 8th-order
            %     y = filtfilt(b,a,y); % Linear-phase implementation
            % 
            %     y   = setdbspl(y,dB_SPL_noise);
            %     innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end
            % 
            %     Wavwrite(innoise,fs,filename);
            % end

            % Generating the test tones

            % Common stim params:
            onset   = t_silence_bef + 115e-3; % temporarily onset equal to the one of the noise + 115 ms
            f       = 1000;
            win     = 1; % 1 = Hanning window

            % Stim 1
            dur     = 5e-3;
            [instim1, t1] = Create_sin4this_exp(f,test_phases(1)*pi,dur,fs,win,onset,dB_SPL,t_total_duration);
            Wavwrite(instim1,fs,filename1);

            for i = 2:N_conditions
                exp1 = sprintf('[instim%.0f, t%.0f] = Create_sin4this_exp(f,%.4f*pi,dur,fs,win,onset,dB_SPL,t_total_duration);',i,i,test_phases(i));
                exp2 = sprintf('Wavwrite(instim%.0f,fs,filename%.0f);',i,i);
                eval(exp1);
                eval(exp2);
            end
            
            nPlots = N_conditions + mod(N_conditions,2);
            figure;
            for i = 1:N_conditions
                subplot( nPlots/2,2,i)
                exp1 = sprintf('plot(t%.0f,instim%.0f); grid on',i,i);
                eval(exp1)
                title(sprintf('initial phase=%.1f * pi',test_phases(i)))
                xlim([onset-dur/4 onset+dur+dur/4]);
            end
            xlabel('t [s]')
            h = Figure2paperfigure(gcf,nPlots/2,2);
            Saveas(h,[Get_TUe_paths('outputs') 'signals-' num2str(options.dB_SPL) '-dB-expII2']);
            
        end

        opts.bPlot      = bPlot;
        
        for i = 1:N_conditions
            
            exp1 = sprintf('out_stim%.0f = Dau1996compare(innoise,instim%.0f,fs,opts);',i,i);
            eval(exp1);
            
            if N_conditions > 3 % remove fields to avoid 'out of memory'
                
                exp3 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig1'');',i,i);
                exp4 = sprintf('out_stim%.0f = Remove_field(out_stim%.0f,''outsig2'');',i,i);
                eval(exp3);
                eval(exp4);
                
            end
            
            exp5 = sprintf('outs.out_stim%.0f = out_stim%.0f;',i,i);
            exp6 = sprintf('clear out_stim%.0f;',i);
            eval(exp5);
            eval(exp6);
            
        end
        
        
    
    case 3 
    %% 3. Signal integration
        % Stimuli:
        %   Noise, white noise, 600 ms, BPF between 20-5000 Hz, 77 dB SPL
        %   Tone 1, 10-ms, 3 kHz, hanning windowed
        %   Tone 2, 20-ms, 3 kHz, hanning windowed
        %   Tone 3, 40-ms, 3 kHz, hanning windowed
        %   The onset of the tones was always 100 ms after masker onset
        options = Ensure_field(options, 'stim_durations',[10 20 40]); % ms
        opts.fc_idx = 3000;
        
        stim_durations      = options.stim_durations;
        N_stim              = length(stim_durations);
        
        infilename      = 'dau1996b_expII3_noisemasker'; % used in Exp 3
        infilename1        = ['dau1996b_expII3_stim01-' num2str(options.dB_SPL)];
        infilename2        = ['dau1996b_expII3_stim02-' num2str(options.dB_SPL)];
        infilename3        = ['dau1996b_expII3_stim03-' num2str(options.dB_SPL)];

        if N_stim > 3
            infilename4        = ['dau1996b_expII3_stim04-' num2str(options.dB_SPL)];
            infilename5        = ['dau1996b_expII3_stim05-' num2str(options.dB_SPL)];
            infilename6        = ['dau1996b_expII3_stim06-' num2str(options.dB_SPL)];
        end

        try 
            [innoise fs] = Wavread([paths.outputs infilename  '.wav']);
        catch
            options.bGenerate   = 1;
            options.bSave_noise = 1;
        end
        
        try
            [instim1 fs] = Wavread([paths.outputs infilename1 '.wav']);
            [instim2 fs] = Wavread([paths.outputs infilename2 '.wav']);
            [instim3 fs] = Wavread([paths.outputs infilename3 '.wav']);

            if N_stim > 3
                [instim4 fs] = Wavread([paths.outputs infilename4 '.wav']);
                [instim5 fs] = Wavread([paths.outputs infilename5 '.wav']);
                [instim6 fs] = Wavread([paths.outputs infilename6 '.wav']);
            end
            options.bGenerate   = 0;
        catch
            options.bGenerate   = 1;
            fs          = 48000;
            options     = Ensure_field(options,'bSave_noise',0);
        end

        options.fs = fs;
        options.typeplot = 2; % Linear scaled

        filename    = [paths.outputs infilename];
        filename1   = [paths.outputs infilename1];
        filename2   = [paths.outputs infilename2];
        filename3   = [paths.outputs infilename3];

        if N_stim > 3
            filename4   = [paths.outputs infilename4];
            filename5   = [paths.outputs infilename5];
            filename6   = [paths.outputs infilename6];
        end

        if options.bGenerate
            t_silence_bef   = 0e-3;
            t_duration      = 600e-3;
            t_silence_aft   = 0e-3;
            t_total_duration = t_silence_bef + t_duration + t_silence_aft;

            Nsil_bef    = round(options.fs*t_silence_bef);
            Nnoise      = round(options.fs*t_duration);
            Nsil_aft    = round(options.fs*t_silence_aft);

            %% Generating the noise

            if options.bSave_noise == 1
                title1 = 'White noise';
                % ymin = -0.15; ymax =  0.15;  yminMU = -100; ymaxMU = 1500;

                % Gen1: white noise, band-pass filtered
                y = wgn(Nnoise,1,1);

                y   =  y(:); % ensures it is a column vector

                Wn = [20 5000]/(options.fs/2); % Normalised cutoff frequency        
                [b,a] = butter(4,Wn); % 8th-order
                y = filtfilt(b,a,y); % Linear-phase implementation

                y   = setdbspl(y,dB_SPL_noise);
                innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end

                Wavwrite(innoise,fs,filename);
            end
            
            % Generating the test tones

            % Common stim params
            onset   = 100e-3;
            f       = 3000;
            win     = 1; % 1 = Hanning window

            % Stim 1
            dur     = 10e-3;
            [instim1, t1] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            Wavwrite(instim1,fs,filename1);

            % Stim 2
            dur     = 20e-3;
            [instim2, t2] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            Wavwrite(instim2,fs,filename2);

            % Stim 3
            dur     = 40e-3;
            [instim3, t3] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
            Wavwrite(instim3,fs,filename3);

            if N_stim > 3
                dur     = stim_durations(4)*1e-3;
                [instim4, t4] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
                Wavwrite(instim4,fs,filename4);

                dur     = stim_durations(5)*1e-3;
                [instim5, t5] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
                Wavwrite(instim5,fs,filename5);

                dur     = stim_durations(6)*1e-3;
                [instim6, t6] = Create_sin4this_exp(f,0,dur,fs,win,onset,dB_SPL,t_total_duration);
                Wavwrite(instim6,fs,filename6);
            end

        end

        opts.bPlot= bPlot; 
        out_stim1 = Dau1996compare(innoise,instim1,fs,opts);
        out_stim2 = Dau1996compare(innoise,instim2,fs,opts);
        out_stim3 = Dau1996compare(innoise,instim3,fs,opts);

        if N_stim > 3 
            out_stim4 = Dau1996compare(innoise,instim4,fs,opts);
            out_stim5 = Dau1996compare(innoise,instim5,fs,opts);
            out_stim6 = Dau1996compare(innoise,instim6,fs,opts);
        end

        %%

        titlefigure = [];

        idx = out_stim1.idx;

        if options.bPlot

            strtitle = sprintf('Noise of %.1f dB, test of %.1f dB\n',options.dB_SPL_noise,options.dB_SPL);
            strtitle1 = [strtitle '(a) signal duration: 10 ms'];
            strtitle2 = [         '(b) signal duration: 20 ms'];
            strtitle3 = [         '(c) signal duration: 40 ms'];
            % Normalised template plot:
            ha = [];
            figure;
            subplot(3,1,1)
            plot(out_stim1.t,out_stim1.template(:,idx)), grid on
            ha(end+1) = gca;
            title(strtitle1)

            subplot(3,1,2)
            plot(out_stim2.t,out_stim2.template(:,idx)), grid on
            ha(end+1) = gca;
            ylabel('normalised amplitude')
            title(strtitle2)

            subplot(3,1,3)
            plot(out_stim3.t,out_stim3.template(:,idx)), grid on
            ha(end+1) = gca;
            xlabel('Time relative to masker onset [s]')
            title(strtitle3)

            h(end+1) = gcf;
            titlefigure{end+1} = sprintf('dau1996b-fig4-sim-noise-%.0fdB-test-%0.fdB',options.dB_SPL_noise,options.dB_SPL);

            AmpMax = 0.06; 

            linkaxes(ha,'xy');
            axis([0 max(out_stim1.t) AmpMax*[-0.25 1]])

            h(end+1) = Figure2paperfigure(h(end),3,1); % new plot in new handle
            close( h(end-1) );
            h(end-1) = [];

            % No-normalised template plot:
            ha = [];
            figure;
            subplot(3,1,1)
            plot(out_stim1.t,out_stim1.template_no_norm(:,idx)), grid on
            ha(end+1) = gca;
            title(strtitle1)

            subplot(3,1,2)
            plot(out_stim2.t,out_stim2.template_no_norm(:,idx)), grid on
            ha(end+1) = gca;
            ylabel('normalised amplitude')
            title(strtitle2)

            subplot(3,1,3)
            plot(out_stim3.t,out_stim3.template_no_norm(:,idx)), grid on
            ha(end+1) = gca;
            xlabel('Time relative to masker onset [s]')
            title(strtitle3)

            h(end+1) = gcf;
            titlefigure{end+1} = sprintf('dau1996b-fig4-sim-no-norm-noise-%.0fdB-test-%0.fdB',options.dB_SPL_noise,options.dB_SPL);

            AmpMax = 250; 

            linkaxes(ha,'xy');
            axis([0 max(out_stim1.t) AmpMax*[-0.25 1]])

            h(end+1) = Figure2paperfigure(h(end),3,1); % new plot in new handle
            close( h(end-1) );
            h(end-1) = [];

            % 
            tp1 = out_stim1.template(:,idx);
            tp2 = out_stim2.template(:,idx);
            tp3 = out_stim3.template(:,idx);
        end

        if N_stim > 3 % remove fields to avoid 'out of memory'
            out_stim1 = Remove_field(out_stim1,'outsig1');
            out_stim1 = Remove_field(out_stim1,'outsig2');
            out_stim2 = Remove_field(out_stim2,'outsig1');
            out_stim2 = Remove_field(out_stim2,'outsig2');
            out_stim3 = Remove_field(out_stim3,'outsig1');
            out_stim3 = Remove_field(out_stim3,'outsig2');
            out_stim4 = Remove_field(out_stim4,'outsig1');
            out_stim4 = Remove_field(out_stim4,'outsig2');
            out_stim5 = Remove_field(out_stim5,'outsig1');
            out_stim5 = Remove_field(out_stim5,'outsig2');
            out_stim6 = Remove_field(out_stim6,'outsig1');
            out_stim6 = Remove_field(out_stim6,'outsig2'); 
        end

        outs.out_stim1 = out_stim1;
        outs.out_stim2 = out_stim2;
        outs.out_stim3 = out_stim3;

        if N_stim > 3
            outs.out_stim4 = out_stim4;
            outs.out_stim5 = out_stim5;
            outs.out_stim6 = out_stim6;
        end

        if options.bSave
            for i=1:length(h)
                try
                    Saveas(h(i), [paths.outputs titlefigure{i}]);
                catch
                    Saveas(h(i), [filename '-handle-' num2str(i)]);
                end
            end
        end
    case 11
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [y,t] = Create_sin4this_exp(f,start_phase,dur,fs,win,onset,SPL,total_duration)

    [y, t]= Create_sin_phase(f,start_phase,dur,fs,win);
    
    % (f,dur,fs,win,onset,dB_SPL_above_thr);
    y  = setdbspl(y,SPL);
    y  = [Gen_silence(onset,fs); y]; 
    
    try % Append silence only if total_duration has been specified
    	y = [y; Gen_silence(total_duration-max(t)-onset-1/fs,fs)];
    end

    t = (1:length(y))/fs; % redefine t
    
end

function t_total_duration = Create_noise(nTag,filename,options)
    
switch nTag
    case 1 | 2 | 10
        
        dB_SPL_noise    = options.dB_SPL_noise;
        fs              = options.fs;
        [t_silence_bef, t_duration, t_silence_after, t_total_duration] = Create_noise_default(nTag);
        
        Nnoise      = round(fs*t_duration);
        
        y = wgn(Nnoise,1,1);

        y   =  y(:); % ensures it is a column vector

        Wn = [20 5000]/(fs/2); % Normalised cutoff frequency        
        [b,a] = butter(4,Wn); % 8th-order
        y = filtfilt(b,a,y); % Linear-phase implementation

        y   = setdbspl(y,dB_SPL_noise);
        innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end

        Wavwrite(innoise,fs,filename);
end
    
end

function [onset, dur, t_silence_aft, t_total_duration] = Create_noise_default(nTag)
    
switch nTag
    case 1 | 2 | 10
        
        onset   = 100e-3;
        dur     = 300e-3;
        t_silence_aft = 200e-3;

        t_total_duration = onset + dur + t_silence_aft;

        Wn = [20 5000]/(fs/2); % Normalised cutoff frequency        
end
    
end