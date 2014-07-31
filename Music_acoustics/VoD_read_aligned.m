function [y_measured y_modelled misc] = VoD_read_aligned(bHPF,info)
% function [y_measured y_modelled misc] = VoD_read_aligned(bHPF,info)
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
%       VoD_read_aligned;
% 
%   % Example 2:
%       bHPF = 1;
%       info.bPlot = 1;
%       [y_measured y_modelled misc] = VoD_read_aligned(bHPF,info);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: VoD_run
% Created on    : 01/07/2014
% Last update on: 31/07/2014 % Update this date manually
% Last use on   : 31/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default inputs:
if nargin == 0
    bHPF = 1;
end

if nargin < 2
    info = [];
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc        = Get_VoD_params(1); % input = 0, to not store 'period' results
misc.hFig   = []; % empty figure handles
misc.near_field_filename = [];
T           = misc.Tmodel;
last_take   = misc.last_take;

numPeriodsSegment = 1; 

subdir_db_vod = Get_TUe_subpaths('db_voice_of_dragon');

% dir calibrated predicted (modelled) signals:
dir_calibrated_m = subdir_db_vod.dir_calibrated_m; 
dir_calibrated_p = subdir_db_vod.dir_calibrated_p; 
dir_f0_m         = subdir_db_vod.dir_f0_m;
dir_f0_p         = subdir_db_vod.dir_f0_p;
dir_meas_def     = subdir_db_vod.dir_meas_def;

% Assumes audio files were already generated:
info = Ensure_field(info,'bPlot',0);
dir_predicted = dir_calibrated_p;
    
for mode = info.modes2check
    
    mode_idx = mode-1;
    modus   = num2str(mode_idx); % 1 to 4
    
    take    = last_take(mode_idx);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Audio files:
    % 1. Measured wav files:
    field = '1';
    filename{1} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    field = '2';
    filename{2} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    file_f0_m   = [dir_f0_m         'modus-' modus '_v' num2str(take) '-' field lblFilter '.txt']; % added on 31/07/2014
    filename{3} = [dir_meas_def 'modus ' modus '_v' num2str(take) '-3.wav'];
    
    % 2. Modelled wav files:
    field = '1';
    filename{4} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
    field = '2';
    filename{5} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
    file_f0_p   = [dir_f0_p         'modus-' modus '-v_' field lblFilter '.txt']; % added on 31/07/2014
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading audio files...
    
    [yfar  fs]  = Wavread(filename{1});
    ynear       = Wavread(filename{2});
    y4period    = Wavread(filename{3});

    [yfarp  fs] = Wavread(filename{4});
    [ynearp fs] = Wavread(filename{5});
    
    misc.near_field_filename{end+1,1} = filename{2}; % measured file name
    misc.near_field_filename{end  ,2} = filename{5}; % modelled file name
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Actual alignment, Version 1: see script version before 31/07/2014
    
    % Actual alignment, Version 2:
    % Truncation according to stored intial time values:
    [ynear , tm] = Do_alignment(  misc.t, ynear , misc.ti_measured(mode_idx)  );
    [yfar  , tm] = Do_alignment(  misc.t, yfar  , misc.ti_measured(mode_idx)  );
    
    [ynearp, tp] = Do_alignment(  misc.t, ynearp, misc.ti_model(mode_idx)     );
    [yfarp , tp] = Do_alignment(  misc.t, yfarp , misc.ti_model(mode_idx)     );
    
    L = min( length(ynear), length(ynearp));
    
    ynear   = Do_truncate(ynear , L);
    yfar    = Do_truncate(yfar  , L);
    ynearp  = Do_truncate(ynearp, L);
    yfarp   = Do_truncate(yfarp , L);
	t       = Do_truncate(tm,L); % same t for measured and modelled

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if info.bPlot
        
        try
            
            n = 4; 
            
            [t_f0m f0m] = Get_F0_praat_from_txt(file_f0_m);
            [f0m tfm] = Do_alignment( t_f0m, f0m, misc.ti_measured(mode_idx) );
            
            [t_f0p f0p] = Get_F0_praat_from_txt(file_f0_p);
            [f0p tfp] = Do_alignment( t_f0p, f0p, misc.ti_model(mode_idx) );
            
            Lf0 = min(length(tfm),length(tfp));
            
            tf = Do_truncate(tfm,Lf0); % same time for both
            
            f0m = Do_truncate(f0m,Lf0);
            f0p = Do_truncate(f0p,Lf0);
            
        catch
            n = 2;
            disp([mfilename '.m: no f0 information found']);
        end
        figure; 
        subplot(n,1,1)
        plot(t,ynear, t, yfar), 
        ha = gca;
        
        legend('near-field','far-field')
        title(sprintf('Aligned time signals, ac.mode = %.0f: measured (above), modelled (bottom)',mode))
        ylabel('Amplitude')
        subplot(n,1,2)
        plot(t,ynearp, t, yfarp), 
        ha(end+1) = gca;
        
        if n == 2
            xlabel('Time [s]')
        end
        
        ylabel('Amplitude')
        
        if n==4
            
            subplot(n,1,3)
            plot(   tf,f0m, ...
                    tf,f0p   )
            grid on, hold on;
            ha(end+1) = gca;
            ylabel('Frequency [Hz]')
            legend('f0 meas','f0 modelled');
            
            subplot(n,1,4)
            plot(   tf,(f0p-f0m)/misc.mf(mode_idx)*100   )
            grid on, hold on;
            ha(end+1) = gca;
            xlabel('Time [s]')
            ylabel('\Delta f / f_n [%]')
            
        end
        
        linkaxes(ha,'x')
        xlim([0 3*misc.Tmodel(mode_idx)]);
        
        misc.hFig(end+1) = gcf;
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    size_yseg_col = round(T(mode_idx)*fs);

    numPeriodsSegment = 1;
    
    disp([mfilename '.m: output structs are going to be buffered containing ' num2str(numPeriodsSegment) ' period per column'])
    disp(['Afterwards, some periods are going to be deleted according to misc.meas_n2, 3, 4 and 5 (for measured data)'])
    
    [y_measured_tmp,non_complete_frame] = buffer(ynear,size_yseg_col*numPeriodsSegment,0);
    exp1 = [' y_measured.yn' num2str(mode) '_short_buf = y_measured_tmp(:,misc.meas_n' num2str(mode) ');'];
    exp2 = [' y_measured.yn' num2str(mode) '_short     = y_measured.yn'       num2str(mode) '_short_buf(:);'];
    
    eval(exp1);
    eval(exp2);
    
    N = eval(['length(y_measured.yn' num2str(mode) '_short);']);
    
    exp3 = ['y_modelled.yn' num2str(mode) '_short  = ynearp( 1:N );']; % Modelled signal having same size than measured signal
    eval(exp3)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    step1 = ['[y_measured.y_seg_f' num2str(mode) ',non_complete_frame] = buffer(yfar,size_yseg_col*numPeriodsSegment,0);'];
    eval(step1);
    
    step1 = ['[y_modelled.y_seg_f' num2str(mode) ',non_complete_frame] = buffer(yfarp,size_yseg_col*numPeriodsSegment,0);'];
    eval(step1);
    
    Exp1 = ['y_measured.yf' num2str(mode) ' = yfar;'];
    eval(Exp1);
    Exp1 = ['y_measured.yn' num2str(mode) ' = ynear;'];
    eval(Exp1);
    
    Exp1 = ['y_modelled.yf' num2str(mode) ' = yfarp;'];
    eval(Exp1);
    Exp1 = ['y_modelled.yn' num2str(mode) ' = ynearp;'];
    eval(Exp1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end