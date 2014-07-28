function VoD_comparisons
% function VoD_comparisons
%
% 1. Description:
%       Make sure you have run previously the script VoD_run.m (without 
%       output parameters) in order to generate up-to-date calibrated 
%       wav-files
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       VoD_comparisons;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on : 28/5/2014
% Last update: 01/07/2014 % Update this date manually
% Last used  : 01/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % tested on 25/06/2014
bPlot_timeseries= 0;    % OK 
bDo_STFT        = 0;    % OK
bDoZwicker      = 0;    % OK
bDo_OB          = 1;    % OK, different buffer size than STFT
bDo_Gammatone   = 0;    % OK, different buffer size than STFT

bHPF            = 1; % 1 = to load band-pass filtered audio files

idx         = 2; % 1 = far field, 2 = near-field
type        = {'f','n'};
info.bSave  = 0;
info.bPlot  = 0;
info.wntype = 1; % Hanning
modes2check = 2:5;

nfft      = 4096; % N-point FFT and STFT
K         = nfft/2; % K effective FFT points
wlen_s    = 15e-3;  % STFT
overlap_p = 5;      % STFT

h       = [];
h_STFT  = [];
h_Loudness = [];
m = [];
s = [];

bDiary = 1;

if bDiary
    
    Diary(mfilename)
    
end

close all;

[y_measured y_modelled misc] = VoD_run(bHPF,info);
fs = misc.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Plot time series:
%       - Last used on: xx/xx/2014
%       - Deleted on  : 01/07/2014
if bPlot_timeseries % To generalise
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show_figures_one_by_one(0.5);

disp([mfilename '.m: analysing ' type{idx} ' field']);

for mode = modes2check
    
    mode_idx = mode-1;
    
    Exp1 = ['lim_max = size(y_modelled.y_seg_' type{idx} num2str(mode_idx) ',2);']; % assuming size(Modelled) = size(Measured)
    eval(Exp1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. STFT:
    %       - Last used on: 30/06/2014
    if bDo_STFT
        wlen = round(wlen_s * fs); % in samples
    
        figure;
        subplot(2,1,1)

        Exp1 = ['stft(y_measured.y' type{idx} num2str(mode_idx) ', fs, nfft, wlen, overlap_p, info.wntype);'];
        eval(Exp1);
        halocal = gca;
        
        title(sprintf('Acoustic mode = %.0f, Measured (top panel), Modelled (bottom panel)',mode))
        xlabel('') % to delete xlabel subplot 2,1,1
        
        cbYLim = [-120 -35];
        hcb = colorbar;
        set(hcb,'YLim',cbYLim);

        subplot(2,1,2)
        Exp1 = ['stft(y_modelled.y' type{idx} num2str(mode_idx) ', fs, nfft, wlen, overlap_p, info.wntype);'];
        eval(Exp1);
        halocal(end+1) = gca;
        hcb = colorbar;
        set(hcb,'YLim',cbYLim);
        
        title('') % to delete title subplot 2,1,2
        
        linkaxes(halocal,'xy');
        
        xlim([0 1])
        ylim([100 1250])
        
        h_STFT(end+1) = gcf;
        
    end
    
    % 3. Zwicker's model
    %       - Last used on: 01/07/2014
    if bDoZwicker
        
        LoudnessI   = figure;
        
        info.bNewFigure = 0;
        info.color = 'b';
        
        % info.bSave = 1;
        info.filename = misc.near_field_filename{mode_idx,1};
        
        Exp1 = ['res{1} = Zwicker_dynamic_loudness_model(y_measured.y' type{idx} num2str(mode_idx) ', fs, info);'];
        eval(Exp1)
        
        info.color = 'r';
        info.filename = misc.near_field_filename{mode_idx,2};
        Exp1 = ['res{2} = Zwicker_dynamic_loudness_model(y_modelled.y' type{idx} num2str(mode_idx) ', fs, info);'];
        eval(Exp1)
        
        subplot(2,1,1)
        legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
        
        subplot(2,1,2)
        legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
        
        info = rmfield(info,'filename');
        % info.bSave = 0;
        
        h_Loudness(end+1) = gcf;
        
        figure;
        ha = subplot(2,1,1); 
        mesh( res{1}.plot_time,res{1}.plot_barkAxis, res{1}.InstantaneousSpecificLoudness )
        
        ha(end+1) = subplot(2,1,2); 
        mesh( res{2}.plot_time,res{2}.plot_barkAxis, res{2}.InstantaneousSpecificLoudness )
        
        Linkaxes_3D(ha);
        xlim([0 20])
        ylim([0 25])
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. OB analysis / one-third OB analysis
    %       - Last used on: 01/07/2014
    if bDo_OB
        
        % Global average for near-field audio files: %%%%%%%%%%%%%%%%%%%%%%
        disp(fprintf('OB-analysis for ac. mode %.0f\n', mode))
        exp1 = ['[P_n' num2str(mode_idx) '_avg_meas, Fc] = Filterbank_analysis( y_measured.yn' num2str(mode_idx) ',fs, 0);'];
        eval(exp1);
        
        exp1 = ['[P_n' num2str(mode_idx) '_avg_meas]'];
        eval(exp1);
        disp('-----')
        exp1 = ['rmsdb( y_measured.yn' num2str(mode_idx) ')'];
        eval(exp1);
        disp('-----')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        info.fs = fs;
        
        figure;
        K = 4096;
        info.bNewFigure = 0;
        
        subplot(3,1,1)
        exp1 = ['y2plot = [y_measured.y' type{idx} num2str(mode_idx) ' y_modelled.y' type{idx} num2str(mode_idx) '];'];
        disp(exp1);
        eval(exp1);
        [N,M] = size(y2plot);
        [win info.wtype] = Get_window(info.wntype,N,M);
        y2plot = y2plot.*win;
        
        if info.bPlot
            
            freqfft(y2plot,K,info);
            try
                legend('Measured','Modelled')
            end
        
            subplot(3,1,2)
            exp1 = ['freqfft(y_measured.y_seg_' type{idx} num2str(mode_idx) ',K,info);'];
            eval(exp1);
            idx2plot = 1:lim_max;
            % idx2plot = 1:4;
            subplot(3,1,3)
            Exp1 = ['freqfft(y_modelled.y_seg_' type{idx} num2str(mode_idx) '(:,idx2plot),K,info);'];
            eval(Exp1);
            Exp1 = ['Print_text_on_plot(num2str( rmsdb(y_modelled.y_seg_' type{idx} num2str(mode_idx) '(:,idx2plot)) ),50);'];
            eval(Exp1);

            h(end+1) = gcf;
            set(h(end),'Position', [1 41 1366 658]);
            
        end
        
        disp('Processing measured file');
        exp1 = ['[Py_' type{idx} num2str(mode_idx) ', Fc] = Filterbank_analysis( y_measured.y_seg_' type{idx} num2str(mode_idx) ',fs, 0);'];
        eval(exp1);
        
        [a b] = eval(['size(Py_' type{idx} num2str(mode_idx) ')']);
        exp1 = ['[Py_' type{idx} num2str(mode_idx) '-repmat(P_n' num2str(mode_idx) '_avg_meas, 1,b)]'];
        
        exp1 = ['Delta_ref_avg = sum( [Py_' type{idx} num2str(mode_idx) '-repmat(P_n' num2str(mode_idx) '_avg_meas, 1,b)] );'];
        eval(exp1);
        
        counter = find(abs(Delta_ref_avg) < 1);
        disp('Cycles that differ in less than 1 dB respect to global band average:')
        disp(['Cycles Nr: ',num2str(counter) ' out of ' num2str(length(Delta_ref_avg)) ' periods'])
        
        disp('Processing modelled file');
        exp1 = ['[Ppy_' type{idx} num2str(mode_idx) ', Fc] = Filterbank_analysis( y_modelled.y_seg_' type{idx} num2str(mode_idx) ',fs, 0);'];
        eval(exp1);
        
        exp1 = ['Delta = Ppy_' type{idx} num2str(mode_idx) ' - Py_' type{idx} num2str(mode_idx) ';'];
        eval(exp1);

        m = [m, transpose( mean( transpose(Delta) ) )];
        s = [m, transpose(  std( transpose(Delta) ) )];

        exp1 = ['Delta' num2str(mode_idx) ' = Ppy_' type{idx} num2str(mode_idx) '- Py_' type{idx} num2str(mode_idx) ';'];
        eval(exp1);

        exp1 = ['[mean2plot,std2plot] = barweb_prepare_data(Py_' type{idx} num2str(mode_idx) ',Ppy_' type{idx} num2str(mode_idx) ');'];
        eval(exp1);

        figure;
        Colors = [1 1 1;0.75 0.75 0.75];
        stPlot.Title = ['Field: ' type{idx} ', Acoustic mode ' num2str(mode_idx+1) ': dB RMS per frequency band'];
        stPlot.XLabel = 'f_c [Hz]'; 
        stPlot.YLabel = 'Relative amplitude [dB]';
        stPlot.SeriesLabel = {'Measured','Modelled'};

        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        ha = gca;
        Set_XTick(ha,round(Fc));
        grid on

        ylabel('RMS value per band [dB]')
        xlabel('Central frequency f_c [Hz]')

        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('Processing measured file');
        exp1 = ['[Pty_n' num2str(mode_idx) ', Fc] = Filterbank_analysis( y_measured.y_seg_n' num2str(mode_idx) ',fs, 1);'];
        eval(exp1);

        disp('Processing modelled file');
        exp1 = ['Ptpy_n' num2str(mode_idx) ' = Filterbank_analysis( y_modelled.y_seg_n' num2str(mode_idx) ',fs, 1);'];
        eval(exp1);

        exp1 = ['[mean2plot,std2plot] = barweb_prepare_data(Pty_n' num2str(mode_idx) ',Ptpy_n' num2str(mode_idx) ')'];
        eval(exp1);

        figure;
        Colors = [1 1 1;0.75 0.75 0.75];
        stPlot.Title = ['Mode ' num2str(mode_idx+1) ': dB RMS per frequency band'];
        stPlot.XLabel = 'f_c [Hz]'; 
        stPlot.YLabel = 'Relative amplitude [dB]';
        stPlot.SeriesLabel = {'Measured','Modelled'};

        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        ha = gca;
        grid on

        Set_XTick(ha,round(Fc))

        ylabel('RMS value per band [dB]')
        xlabel('Central frequency f_c [Hz]')

        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);
        % Show_figures_one_by_one(0.5);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. Gammatone
    %       - Last used on: xx/xx/2014
    if bDo_Gammatone

        disp('Processing measured file');
        exp1 = ['[Pgy_n' num2str(mode_idx) ', Fc] = Filterbank_analysis( y_measured.y_seg_n' num2str(mode_idx) ',fs, 2);'];
        eval(exp1);

        disp('Processing modelled file');
        exp1 = ['[Pgpy_n' num2str(mode_idx) ', Fc] = Filterbank_analysis( y_modelled.y_seg_n' num2str(mode_idx) ',fs, 2);'];
        eval(exp1);

        exp1 = ['[mean2plot,std2plot] = barweb_prepare_data(Pgy_n' num2str(mode_idx) ',Pgpy_n' num2str(mode_idx) ');'];
        eval(exp1);

        figure;
        Colors = [1 1 1;0.75 0.75 0.75];
        stPlot.Title = ['Mode ' num2str(mode_idx+1) ': dB RMS per frequency band (Gammatone)'];
        stPlot.XLabel = 'f_c [Hz]'; 
        stPlot.YLabel = 'Relative amplitude [dB]';
        stPlot.SeriesLabel = {'Measured','Modelled'};

        barweb(mean2plot,std2plot,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
        ha = gca;
        Set_XTick(ha,round(Fc));
        grid on

        ylabel('RMS value per band [dB]')
        xlabel('Central frequency f_c [Hz]')

        h(end+1) = gcf;
        set(h(end),'Position', [1 41 1366 658]);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info.bSave
    
    paths.outputs   = Get_TUe_paths('outputs');
    
    if isfield(misc,'hFig')
        for j = 1:length(misc.hFig)
            Saveas(misc.hFig(j),[paths.outputs 'aligned-ac_mode-' num2str(j+1)]);
        end
    end
    
    if bDo_STFT
        for j = 1:length(h_STFT)
            Saveas(h_STFT(j),[paths.outputs 'stft-ac_mode-' num2str(j+1)]);
        end
    end
    
    for j = 1:length(h)
        Saveas(h(j)     ,[paths.outputs 'freq_analysis-' num2str(j)]);
    end
    
end

if bDiary
    diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])