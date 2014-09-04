function [misc] = VoD_write_aligned(info)
% function [misc] = VoD_write_aligned(info)
%
% 1. Description:
%       Similar to VoD_run
% 
%       Calibration: measured near-field file is the reference, dBSPL [dB]. 
%       Same level is used for near-field predicted audio file. Delta dB is
%       determined for near to far-field wav files and the same attenuation
%       is applied to predicted far-field file
%       
%       In terms of alignment, predicted audio files will be the reference
% 
%       Outputs are structs with far and near fields per modus 1-4 (modes 2 
%       to 5):  y_measured 
%               y_modelled
%       Additionally:
%               misc.fs     - Audio sampling frequency
%               misc.Tmodel - rotation period for VoD model according to Hirschberg2013
%               misc.near_field_filename - wav file names read. Column 1 is
%                             for measured file names, column 2 for modelled 
%                             file names.
% 2. Additional info:
%   Tested cross-platform: Yes
%
% 3. Stand-alone example:
%   % Example 1: to get all the plots, for measured and predicted VoD files
%       VoD_write_aligned;
% 
%   % Example 2:
%       info.Dest_folder = Get_TUe_paths('outputs');
%       info.bPlot = 1;
%       VoD_read_aligned(info);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: VoD_read_aligned
% Created on    : 29/08/2014
% Last update on: 29/08/2014 % Update this date manually
% Last use on   : 29/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default inputs:
bHPF = 1;

if nargin == 0
    info = [];
end

info = Ensure_field(info,'dest_folder',Get_TUe_paths('outputs'));
dest_folder = info.dest_folder;

if bHPF
    disp([mfilename '.m: HPF is going to be applied (calibration takes this into account, just be aware)'])
    pause(2)
end

if bHPF
    lblFilter = 'filt'; % Suffix for wav-files, in case they were filtered
else
    lblFilter = '';
end

info = Ensure_field(info,'modes2check',2:5);
info = Ensure_field(info,'bSave',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc        = Get_VoD_params(1); % input = 0, to not store 'period' results
misc.hFig   = []; % empty figure handles
misc.ha     = [];
misc.near_field_filename = [];
T           = misc.Tmodel;
last_take   = misc.last_take;

subdir_db_vod = Get_TUe_subpaths('db_voice_of_dragon');

% dir calibrated predicted (modelled) signals:
dir_calibrated_m = subdir_db_vod.dir_calibrated_m; 
dir_calibrated_p = subdir_db_vod.dir_calibrated_p; 
dir_calibrated_ms = subdir_db_vod.dir_calibrated_ms; 
dir_calibrated_ps = subdir_db_vod.dir_calibrated_ps; 
dir_f0_m         = subdir_db_vod.dir_f0_m;
dir_f0_p         = subdir_db_vod.dir_f0_p;
dir_meas_def     = subdir_db_vod.dir_meas_def;
dir_wav_all      = subdir_db_vod.dir_wav_all;

% Assumes audio files were already generated:
info = Ensure_field(info,'bPlot',0);
dir_predicted = dir_calibrated_p;
    
for acmode = info.modes2check
    
    mode_idx = acmode-1;
    modus   = num2str(mode_idx); % 1 to 4
    
    take    = last_take(mode_idx);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Audio files:
    % 1. Measured wav files:
    field = '1';
    filename{1} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    field = '2';
    filename1 = ['modus-' modus '_v' num2str(take) '-' field lblFilter];
    filename{2} = [dir_calibrated_m filename1 '.wav'];
    file_f0_m   = [dir_f0_m         filename1 '.txt']; % added on 31/07/2014
    filename{3} = [dir_meas_def     'modus ' modus '_v' num2str(take) '-3.wav'];
    filename1new = sprintf('meas-ac-mode-%.0f',acmode);
    
    % 2. Modelled wav files:
    field = '1';
    filename{4} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
    
    field = '2';
    filename2 = ['modus-' modus '-v_' field lblFilter];
    filename{5} = [dir_calibrated_p filename2 '.wav'];
    file_f0_p   = [dir_f0_p         filename2 '.txt']; % added on 31/07/2014
    filename2new = sprintf('model-ac-mode-%.0f',acmode);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading audio files...
    
    [ynear, fs] = Wavread(filename{2});
    y4period    = Wavread(filename{3});

    [ynearp fs] = Wavread(filename{5});
    
    misc.near_field_filename{end+1,1} = filename{2}; % measured file name
    misc.near_field_filename{end  ,2} = filename{5}; % modelled file name
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Actual alignment, Version 1: see script version before 31/07/2014
    
    % Actual alignment, Version 2:
    % Truncation according to stored intial time values:
    [ynear , tm] = Do_alignment(  misc.t, ynear , misc.ti_measured(mode_idx)  );
        
    [ynearp, tp] = Do_alignment(  misc.t, ynearp, misc.ti_model(mode_idx)     );
    
    L = min( length(ynear), length(ynearp));
    
    ynear   = Do_truncate(ynear , L);
    ynearp  = Do_truncate(ynearp, L);
    t       = Do_truncate(tm,L); % same t for measured and modelled

    misc.Excerpt_m{mode_idx} = [dest_folder filename1new];
    misc.Excerpt_p{mode_idx} = [dest_folder filename2new]; 
    
    if info.bSave
        info = Ensure_field(info, 'time2save', 5);
        
        disp('Include CALIBRATION step!')
        ynear2 = resample(ynear,44100,fs);
        ynear2 = ynear2(1:round(info.time2save*44100));
        Wavwrite(ynear2,44100,misc.Excerpt_m{mode_idx});
        Get_F0_AC_praat([misc.Excerpt_m{mode_idx} '.wav'],[misc.Excerpt_m{mode_idx} '.txt']);
        
        disp('Include CALIBRATION step!')
        ynearp2 = resample(ynearp,44100,fs);
        ynearp2 = ynearp2(1:round(info.time2save*44100));
        Wavwrite(ynearp2,44100,misc.Excerpt_p{mode_idx});
        Get_F0_AC_praat([misc.Excerpt_p{mode_idx} '.wav'],[misc.Excerpt_p{mode_idx} '.txt']);
        
        disp([mfilename '.m: Synchronised file names returned...'])
    else
        disp([mfilename '.m: Excerpt file names returned but not saved'])
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if info.bPlot
        
        try
            
            n = 4; 
            
            [tfm f0m] = Get_F0_praat_from_txt( [misc.Excerpt_m{mode_idx} '.txt'] );
            
            [tfp f0p] = Get_F0_praat_from_txt( [misc.Excerpt_p{mode_idx} '.txt'] );
            
            Lf0 = min(length(tfm),length(tfp));
            
            tf = Do_truncate(tfm,Lf0); % same time for both
            
            f0m = Do_truncate(f0m,Lf0);
            f0p = Do_truncate(f0p,Lf0);
            
            Exp1 = sprintf('misc.tf0%.0f = tf' ,mode_idx);
            Exp2 = sprintf('misc.f0p%.0f = f0p;',mode_idx);
            Exp3 = sprintf('misc.f0m%.0f = f0m;',mode_idx);
            eval(Exp1)
            eval(Exp2)
            eval(Exp3)
        catch
            n = 2;
            disp([mfilename '.m: no f0 information found']);
            disp('set info.bSave to 1 and re-run this script')
            pause(2)
        end
        
        figure; 
        subplot(n,1,1)
        plot(t,ynear)
        ha = gca;
        
        title(sprintf('Aligned time signals, ac.mode = %.0f: measured (above), modelled (bottom)',acmode))
        ylabel('Amplitude')
        subplot(n,1,2)
    
        plot(t,ynearp,'r'), 
        ha(end+1) = gca;
        
        if n == 2
            xlabel('Time [s]')
        end
        
        ylabel('Amplitude')
        
        if n==4
            
            subplot(n,1,3)
            plot(   tf,f0m, 'b-','LineWidth',1), grid on, hold on
            plot(   tf,f0p, 'r-','LineWidth',2)
            grid on, hold on;
            ha(end+1) = gca;
            ylabel('Frequency [Hz]')
            title('Estimated fundamental frequency f0')
                        
            subplot(n,1,4)
            plot(   tf,(f0p-f0m)/misc.mf(mode_idx)*100   ,'k')
            grid on, hold on;
            title('\Delta f0')
            ha(end+1) = gca;
            xlabel('Time [s]')
            ylabel('\Delta f / f_n [%]')
            
        end
        
        linkaxes(ha,'x')
        xlim([0 3*misc.Tmodel(mode_idx)]);
        
        misc.hFig(end+1) = gcf;
        misc.ha(end+1) = ha(end);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Commented on 3/09/2014, the following lines were generating the 
        % modified plot version normalised to JND, see figure 14, report on 29/08/2014
        % try 
            % figure;
            % misc.hFig(end+1) = gcf;
            % CloneFig(misc.hFig(end-1),misc.hFig(end));
            % 
            % subplot(n,1,4)
            % JND = Get_JND_freq_tones(misc.mf(mode_idx));
            % f_JND = misc.mf(mode_idx) + JND/2;
            % xlim([0 3*misc.Tmodel(mode_idx)]);
            % 
            % plot(   [min(tf) max(tf)],[ 100  100],'g--'), hold on
            % Diff_f = (f0p-f0m)/JND*100;
            % plot(   tf, abs(Diff_f),'k')
            % Diff_f = abs(Delete_NaN_columns(Diff_f));
            % p05 = percentile(Diff_f,100- 5);
            % p10 = percentile(Diff_f,100-10);
            % p20 = percentile(Diff_f,100-20);
            % p50 = percentile(Diff_f,100-50);
            % p90 = percentile(Diff_f,100-90);
            % p95 = percentile(Diff_f,100-95);
            % 
            % grid on;
            % title(['\Delta f0, normalised to JND [%]. (Percentiles [%]: P_{5,10,90,95} =' sprintf(' %.0f, %.0f, %.1f, %.1f)',p05,p10,p90,p95)]);
            % ha(end+1) = gca;
            % xlabel('Time [s]')
            % ylabel('norm( \Delta f ) [%]')
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end