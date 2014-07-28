function [y_measured y_modelled misc] = VoD_run(bHPF,info,take,modes2check)
% function [y_measured y_modelled misc] = VoD_run(bHPF,info,take,modes2check)
%
% 1. Description:
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
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   % Example 1: to get all the plots, for measured and predicted VoD files
%   %            It does generate the audio files
% 
%   VoD_run;
% 
%   % Example 2: to get plots, for measured and predicted VoD files, mode = 2
%   mode = 2;
%   VoD_run(0,[],mode);
% 
%   % Example 3: to get plots, for measured (take 3) and predicted VoD files, mode = 4
%   mode = 4;
%   take = 3;
%   VoD_run(0,[],mode,take);
%  
%   % Example 4:
%    bHPF = 1;
%    [y_measured y_modelled misc] = VoD_run(bHPF);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/05/2014
% Last update on: 27/07/2014 % Update this date manually
% Last use on   : 27/07/2014 % Update this date manually
% 
% Original file name: Run_voice_of_dragon, changed on: 27/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default inputs:
if nargin == 0
    bHPF = 1;
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

if nargin < 4
    modes2check = 2:5; % 2 to 5
end

if nargin < 2
    info = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc        = Get_VoD_params(0); % input = 0, to not store 'period' results
misc.hFig   = []; % empty figure handles
misc.near_field_filename = [];
T           = misc.Tmodel;
last_take   = misc.last_take;

numPeriodsSegment = 1; 
% y_seg   = [];
% yp_seg  = [];

subdir_db_vod = Get_TUe_subpaths('db_voice_of_dragon');

% dir calibrated predicted (modelled) signals:
dir_calibrated_m = subdir_db_vod.dir_calibrated_m; 
dir_calibrated_p = subdir_db_vod.dir_calibrated_p; 
dir_meas_def     = subdir_db_vod.dir_meas_def;

% info.bSave = input('Type 1 to save calibrated audio files. Type 0 otherwise: ');
if nargout ~=0 % Assumes audio files were already generated
    
    info = Ensure_field(info,'bSave',0); 
    info = Ensure_field(info,'bPlot',0);
    dir_predicted = dir_calibrated_p;
    
else
    
    Mkdir(dir_calibrated_m);
    Mkdir(dir_calibrated_p);
    
    info.bSave      = input('Do you want to save your results?: ');
    info            = Ensure_field(info,'bPlot',0);
    dir_predicted   = subdir_db_vod.dir_predicted_txt;
    
    disp(['Modelled audio files being taken from: ' dir_predicted])
    disp(['Is this OK?...'])
    disp('Press any button to continue...')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause;
    
end

info.bRead = ~info.bSave;

disp(['Reading measured files from: ' dir_meas_def ])
disp(['Reading modelled files from: ' dir_predicted])

for mode = modes2check
    
    mode_idx = mode-1;
    modus   = num2str(mode_idx); % 1 to 4
    
    if nargin <= 2 
        take    = last_take(mode_idx);
    elseif length(take) == 0
        take    = last_take(mode_idx);
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Audio files:
    % 1. Measured wav files:
    if info.bRead == 0 % then wav-files do not exist yet
        field = '1';
        filename{1} = [dir_meas_def     'modus ' modus '_v' num2str(take) '-' field '.wav'];
        field = '2';
        filename{2} = [dir_meas_def     'modus ' modus '_v' num2str(take) '-' field '.wav'];
    else
        field = '1';
        filename{1} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
        field = '2';
        filename{2} = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
    end
    filename{3} = [dir_meas_def 'modus ' modus '_v' num2str(take) '-3.wav'];
    
    % 2. Modelled wav files:
    if info.bRead == 0 % then wav-files do not exist yet
        field = '1';
        filename{4} = [dir_predicted 'modus-' modus '-v_' field '.wav'];
        field = '2';
        filename{5} = [dir_predicted 'modus-' modus '-v_' field '.wav'];
    else
        field = '1';
        filename{4} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
        field = '2';
        filename{5} = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
    end
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
    % Truncation according to stored intial time values:
    t       = misc.t;
    idx_t  = min(find( t >= misc.ti_measured(mode_idx) ));
    idx_tp = min(find( t >=    misc.ti_model(mode_idx) ));
    
    yneartrunc   = ynear(idx_t:end);
    yfartrunc    = yfar(idx_t:end);
    ynearptrunc  = ynearp(idx_tp:end);
    yfarptrunc   = yfarp(idx_tp:end);
    
    L = min( length(yneartrunc), length(ynearptrunc));
    
    yneartrunc   = Do_truncate(yneartrunc , L);
    yfartrunc    = Do_truncate(yfartrunc  , L);
    ynearptrunc  = Do_truncate(ynearptrunc, L);
    yfarptrunc   = Do_truncate(yfarptrunc , L);
	ttrunc       = Do_truncate(t,L);
    
    if info.bPlot
        
        figure; 
        subplot(2,1,1)
        plot(ttrunc,yneartrunc, ttrunc, yfartrunc), 
        xlim([0 1])
        legend('near-field','far-field')
        title(sprintf('Aligned time signals, ac.mode = %.0f: measured (above), modelled (bottom)',mode))
        ylabel('Amplitude')
        subplot(2,1,2)
        plot(ttrunc,ynearptrunc, ttrunc, yfarptrunc), 
        xlim([0 1])
        xlabel('Time [s]')
        ylabel('Amplitude')

        misc.hFig(end+1) = gcf;
        
    end
    
    if info.bSave == 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HPF with default fc = 100 Hz
        if bHPF == 1
            info.fs = fs;
            fch = 100;
            % fcl = 
            ynear   = freqfftwhpf(ynear ,info,fch); % default fc = 50 Hz;
            ynear   = freqfftwlpf(ynear ,info);
            
            ynearp  = freqfftwhpf(ynearp,info,fch);
            ynearp  = freqfftwlpf(ynearp,info);

            yfar    = freqfftwhpf(yfar  ,info,fch);
            yfar    = freqfftwlpf(yfar  ,info);

            yfarp   = freqfftwhpf(yfarp ,info,fch);
            yfarp   = freqfftwlpf(yfarp ,info);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calibration, only in case of file generation

        dBSPL = dbspl(ynear);
        Att2farfield = dBSPL - dbspl(yfar);

        ynearp = setdbspl(ynearp, dBSPL);
        yfarp  = setdbspl(yfarp , dBSPL-Att2farfield);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Saving data (after truncation)
        N = length(ynear );
        M = length(ynearp);
        
        P = min(N,M); % they are supposed to be the same value
        ynear   = ynear(1:P);
        ynearp  = ynearp(1:P);
        yfar    = yfar(1:P);
        yfarp   = yfarp(1:P);
        
        field = '2';
        outfilename = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
        Wavwrite(ynear,fs,outfilename);
        
        field = '1';
        outfilename = [dir_calibrated_m 'modus-' modus '_v' num2str(take) '-' field lblFilter '.wav'];
        Wavwrite(yfar,fs,outfilename);
        
        field = '2';
        outfilename = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
        Wavwrite(ynearp,fs,outfilename);    

        field = '1';
        outfilename = [dir_calibrated_p 'modus-' modus '-v_' field lblFilter '.wav'];
        Wavwrite(yfarp,fs,outfilename);

    else
        
        disp('Assigning truncated data to data output...')
        ynear   = yneartrunc;
        yfar    = yfartrunc;
        ynearp  = ynearptrunc;
        yfarp   = yfarptrunc;
        %misc.t  = ttrunc;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    size_yseg_col = round(T(mode_idx)*fs);
    
    step1 = ['y_measured.y_seg_n' num2str(mode_idx) ' = buffer(ynear,size_yseg_col*numPeriodsSegment,0);'];
    step2 = ['y_measured.y_seg_n' num2str(mode_idx) '(:,end)=[];'];
    eval(step1);
    eval(step2);

    step1 = ['y_modelled.y_seg_n' num2str(mode_idx) ' = buffer(ynearp,size_yseg_col*numPeriodsSegment,0);'];
    step2 = ['y_modelled.y_seg_n' num2str(mode_idx) '(:,end)=[];'];
    eval(step1);
    eval(step2);

    step1 = ['y_measured.y_seg_f' num2str(mode_idx) ' = buffer(yfar,size_yseg_col*numPeriodsSegment,0);'];
    step2 = ['y_measured.y_seg_f' num2str(mode_idx) '(:,end)=[];'];
    eval(step1);
    eval(step2);

    step1 = ['y_modelled.y_seg_f' num2str(mode_idx) ' = buffer(yfarp,size_yseg_col*numPeriodsSegment,0);'];
    step2 = ['y_modelled.y_seg_f' num2str(mode_idx) '(:,end)=[];'];
    eval(step1);
    eval(step2);
    
    Exp1 = ['y_measured.yf' num2str(mode_idx) ' = yfar;'];
    eval(Exp1);
    Exp1 = ['y_measured.yn' num2str(mode_idx) ' = ynear;'];
    eval(Exp1);
    
    Exp1 = ['y_modelled.yf' num2str(mode_idx) ' = yfarp;'];
    eval(Exp1);
    Exp1 = ['y_modelled.yn' num2str(mode_idx) ' = ynearp;'];
    eval(Exp1);
    
end

if nargin ~= 0 % then we display our new structs...
    
    misc.fs = fs;
    y_measured
    y_modelled
    
    try
        disp('RMS values for far field (meas,modelled) and near-field (meas,modelled), respectively. Ascending mode number');
        disp('mode 2')
        disp([rmsdb(y_measured.yf1) rmsdb(y_modelled.yf1) rmsdb(y_measured.yn1) rmsdb(y_modelled.yn1)])
        disp('mode 3')
        disp([rmsdb(y_measured.yf2) rmsdb(y_modelled.yf2) rmsdb(y_measured.yn2) rmsdb(y_modelled.yn2)])
        disp('mode 4')
        disp([rmsdb(y_measured.yf3) rmsdb(y_modelled.yf3) rmsdb(y_measured.yn3) rmsdb(y_modelled.yn3)])
        disp('mode 5')
        disp([rmsdb(y_measured.yf4) rmsdb(y_modelled.yf4) rmsdb(y_measured.yn4) rmsdb(y_modelled.yn4)])
    catch
        warning('Disregard this warning if you just processed partially the VoD files...');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end