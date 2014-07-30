function VoD_comparisons2(bDiary)
% function VoD_comparisons2(bDiary)
%
% 1. Description:
%       Same than VoD_comparisons, but here analysis is done over aligned/
%       truncated data.
% 
%       Make sure you have run previously the script VoD_run.m (without 
%       output parameters) in order to generate up-to-date calibrated 
%       wav-files
% 
% 2. Additional info:
%       Tested cross-platform: No
%       near-field: YES, OK
%       far-field : NO,  OK 
%
% 3. Stand-alone example:
%       bDiary = 1; % to generate a log-file
%       VoD_comparisons2(bDiary);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 01/07/2014
% Last update on: 28/07/2014 % Update this date manually
% Last use on   : 28/07/2014 % Update this date manually
% 
% Original file name: VoD_comparisons.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    bDiary = 1;
end
Diary(mfilename,bDiary)

%                         % tested on 01/07/2014
bDo_STFT        = 1;    % Analysis 1
bDoZwicker      = 0;    % Analysis 2
bMIR            = 0;    % Analysis discarded according to meeting on 23/07/2014

% bDo_OB          = 1;    % different buffer size than STFT
% bDo_Gammatone   = 0;    % different buffer size than STFT

bHPF        = 1; % 1 = to load band-pass filtered audio files
idx         = 2; % 1 = far field, 2 = near-field
type        = {'f','n'};
info.bSave  = 1;
info.bPlot  = 1;
info.modes2check = 2:5;
bGridOn     = 1;

info.wntype = 1;    % Hanning
nfft      = 4096;   % STFT, N-point FFT and STFT
K         = nfft/2; % STFT, K effective FFT points
wlen_s    = 15e-3;  % STFT
overlap_p = 5;      % STFT

h_STFT      = [];
h_STFT_Mesh = [];
h_Loudness  = [];
% m = [];
% s = [];
 
close all;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Signal alignment:
%       - You get the following plots (if info.bPlot == 1):
%           * 4 x Near and far field aligned signals (1 per acoustic mode)
[ymeasured ymodelled misc] = VoD_read_aligned(bHPF,info);
fs = misc.fs;
 
% Show_figures_one_by_one(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' type{idx} ' field']);
 
for mode = info.modes2check
    
    mode_idx = mode-1;
    ymeas       = eval(['ymeasured.y' type{idx} num2str(mode) '_short;']);
    ymodel      = eval(['ymodelled.y' type{idx} num2str(mode) '_short;']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. STFT:
    %       - Last used on: 28/07/2014
    if bDo_STFT
        
        YLim_spec = [200 700; 400 900; 600 1200; 600 1200]; % manually determined
        info.normalise_time_factor = misc.Tmodel(mode_idx);
        
        wlen = round(wlen_s * fs); % in samples
    
        figure;
        subplot(2,1,1)

        stft(ymeas, fs, nfft, wlen, overlap_p, info.wntype,info);
        [Hmeas, f, t] = stft(ymeas, fs, nfft, wlen, overlap_p, info.wntype);
        HdBmeas = 20*log10(abs(Hmeas));
        
        halocal = gca;
        
        title(sprintf('Acoustic mode = %.0f, Measured (top panel), Modelled (bottom panel)',mode))
        xlabel('') % to delete xlabel subplot 2,1,1
        
        cbYLim = [-120 -35];
        hcb = colorbar;
        set(hcb,'YLim',cbYLim);

        subplot(2,1,2)
        stft(ymodel, fs, nfft, wlen, overlap_p, info.wntype,info);
        Hmodel   = stft(ymodel, fs, nfft, wlen, overlap_p, info.wntype);
        HdBmodel = 20*log10(abs(Hmodel));
        halocal(end+1) = gca;
        hcb = colorbar;
        set(hcb,'YLim',cbYLim);
        
        title('') % to delete title subplot 2,1,2
        
        linkaxes(halocal,'xy');
        
        xlim([0 3*misc.Tmodel(mode_idx)/info.normalise_time_factor]) % showing 3 periods
        ylim([YLim_spec(mode_idx,:)])
        
        h_STFT(end+1) = gcf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        if bGridOn
            
            if mode_idx == 1 % Only done once!
                Delta_coil_ti = nan(4,100);
                
                bGetT = 1;
                figure;
                tmp = Get_VoD_params(bGetT);
                t_coil_noalign = tmp.t_coil; % Correct this with the use of .t_coil
                t_coil = t_coil_noalign-repmat(tmp.ti_measured',1,length(t_coil_noalign));
                t_coil( find(t_coil < 0) ) = NaN;
                t_coil = Delete_NaN_columns(t_coil);
                close;
                close;
            end
            
            CurrT = misc.Tmodel(mode_idx);
            pos     = 0;% - misc.ti_model(mode_idx);
            
            Freq_res = misc.mf(mode_idx);
            
            % Predominant frequencies:
            fpmeas = Get_predominant_values(f,HdBmeas);
            subplot(2,1,1)
            hold on;
            plot(get(gca,'XLim'),[Freq_res Freq_res],'g--')
            plot(t/info.normalise_time_factor,fpmeas,'r')

            subplot(2,1,2)
            hold on;
            fpmodel = Get_predominant_values(f,HdBmodel);
            hold on;
            plot(get(gca,'XLim'),[Freq_res Freq_res],'g--')
            plot(t/info.normalise_time_factor,fpmodel,'r')

            % Cross initial position:
            %   - Assumes perfect synchronisation between measured and modelled
            misc.ti_measured(mode_idx);
            
            Odd     = 1;
            Lim_sup = max(get(gca,'XLim'));
            while pos(end) < Lim_sup % Will draw lines up to upper limit of the axis
                
                if Odd % Real Ti
                    subplot(2,1,1)
                    plot([pos(end) pos(end)], [0 1500],'LineWidth',2)
                    
                    subplot(2,1,2)
                    plot([pos(end) pos(end)], [0 1500],'LineWidth',2)
                else
                    subplot(2,1,1)
                    plot([pos(end) pos(end)], [0 1500],'b--','LineWidth',1)
                    
                    subplot(2,1,2)
                    plot([pos(end) pos(end)], [0 1500],'b--','LineWidth',1)
                end
                Odd = ~Odd;
                pos(end+1) = pos(end) + CurrT/2/info.normalise_time_factor;
                
            end
            
            t_coil_this_mode = Delete_NaN_columns( t_coil(mode_idx,:) );
            t_coil_this_mode = t_coil_this_mode(find( t_coil_this_mode > pos(1)))/info.normalise_time_factor; 
            t_coil_this_mode = t_coil_this_mode(find( t_coil_this_mode < Lim_sup + misc.Tmodel(mode_idx))); 
            
            Delta_tmp = t_coil_this_mode; % - pos(1:2:Lim_sup);
            Delta_coil_ti(mode_idx, 1:length(Delta_tmp)) = Delta_tmp;
            
            for j = 1:length(Delta_tmp)
                subplot(2,1,1)
                plot([Delta_tmp(j) Delta_tmp(j)], [0 1500],'g--','LineWidth',2)
                subplot(2,1,2)
                plot([Delta_tmp(j) Delta_tmp(j)], [0 1500],'g--','LineWidth',2)
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate plots either in 3D or 2D (just stored, not shown):
        tmpOption.bPlot3D = 0;
        
        if tmpOption.bPlot3D
            figure;
            subplot(2,1,1)
        end
        
        Mesh(t/info.normalise_time_factor,f,HdBmeas,tmpOption);
        
        if tmpOption.bPlot3D
            set(gca,'CLim',[min(min(HdBmeas))-40 max(max((HdBmeas)))+20])
            view([110 30])
            colorbar
        end
        
        if tmpOption.bPlot3D
            subplot(2,1,2)
        end
        
        Mesh(t/info.normalise_time_factor,f,HdBmodel,tmpOption);
        
        if tmpOption.bPlot3D
            set(gca,'CLim',[min(min(HdBmeas))-40 max(max((HdBmeas)))+20]) % scaled respect to the other plot
            view([110 30])
            colorbar
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h_STFT_Mesh(end+1) = gcf;
        
    end
    
    % 3. Zwicker's model
    %       - Last used on: 03/07/2014
    if bDoZwicker
        
        figure;
        info.bNewFigure = 0;
        info.color = 'b';
        
        info.filename = misc.near_field_filename{mode_idx,1};
        
        info.bSave = 0;
        res{1} = Zwicker_dynamic_loudness_model(ymeas, fs, info);
        
        info.color = 'r';
        info.filename = misc.near_field_filename{mode_idx,2};
        res{2} = Zwicker_dynamic_loudness_model(ymodel, fs, info);
        info.bSave = 1;
        
        subplot(2,1,1)
        legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
        halocal = gca;
        
        scale = 1.8;
        Ylimit = get(halocal(end),'YLim');
        set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);
        
        subplot(2,1,2)
        legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
        halocal(end+1) = gca;
        scale = 1.8;
        Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
        set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);
        
        linkaxes(halocal,'x')
        
        info = rmfield(info,'filename');
        % info.bSave = 0;
        
        h_Loudness(end+1) = gcf;
        
        if mode_idx == 4

            ii = 1;
            figure;
            info.bNewFigure = 0;
            info.color = 'b';

            Nfinal = round(misc.Tmodel(ii)*2*44100); %  seconds

            filename1 = 'modus-1_v2-2filt-fc-251-Hz.wav';
            [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
            ymeas2 = ymeas2(1:Nfinal);
            info.bSave = 0;
            res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

            filename2 = 'modus-1-v_2filt-fc-251-Hz.wav';
            [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
            ymodel2 = ymodel2(1:Nfinal);
            info.color = 'r';
            res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
            info.bSave = 1;

            subplot(2,1,1)
            legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
            halocal = gca;
            title(name2figname(filename1))

            scale = 1.8;
            Ylimit = get(halocal(end),'YLim');
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

            subplot(2,1,2)
            legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
            halocal(end+1) = gca;
            Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);
            hZwicker_band = gcf;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % High frequency band

            figure;

            filename1 = 'modus-1_v2-2filt-fc-1000-Hz.wav';
            [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
            ymeas2 = ymeas2(1:Nfinal);
            info.bSave = 0;
            info.color = 'b';
            res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

            filename2 = 'modus-1-v_2filt-fc-1000-Hz.wav';
            [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
            ymodel2 = ymodel2(1:Nfinal);
            info.color = 'r';
            res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
            info.bSave = 1;

            subplot(2,1,1)
            legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
            halocal(end+1) = gca;
            title(name2figname(filename1))

            scale = 1.8;
            Ylimit = get(halocal(end),'YLim');
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

            subplot(2,1,2)
            legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
            halocal(end+1) = gca;
            Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

            hZwicker_band(end+1) = gcf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % High frequency band

            figure;

            filename1 = 'modus-1_v2-2filt-fc-3981-Hz.wav';
            [ymeas2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\' filename1]);
            ymeas2 = ymeas2(1:Nfinal);
            info.bSave = 0;
            info.color = 'b';
            res{1} = Zwicker_dynamic_loudness_model(ymeas2, fs, info);

            filename2 = 'modus-1-v_2filt-fc-3981-Hz.wav';
            [ymodel2 fs] = Wavread([Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\' filename2]);
            ymodel2 = ymodel2(1:Nfinal);
            info.color = 'r';
            res{2} = Zwicker_dynamic_loudness_model(ymodel2, fs, info);
            info.bSave = 1;

            subplot(2,1,1)
            legend(['meas. ' res{1}.txt_loudness],['modelled ' res{2}.txt_loudness])
            halocal(end+1) = gca;
            title(name2figname(filename1))

            scale = 1.8;
            Ylimit = get(halocal(end),'YLim');
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

            subplot(2,1,2)
            legend(['meas.' res{1}.txt_db],['modelled ' res{2}.txt_db])
            halocal(end+1) = gca;
            Ylimit = get(halocal(end),'YLim'); %YLim = [-0.04 0.04]
            set(halocal(end),'YLim',[Ylimit(1) Ylimit(2)*scale]);

            linkaxes(halocal,'x')

            hZwicker_band(end+1) = gcf;
        end 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. MIR analysis:
    %       - Analysis discarded according to meeting on 23/07/2014
    if bMIR
        VoD_MIRtoolbox(ymeas , fs);
        VoD_MIRtoolbox(ymodel, fs);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

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
        for j = 1:length(h_STFT_Mesh)
            Saveas(h_STFT_Mesh(j),[paths.outputs 'stft-mesh-ac_mode-' num2str(j+1)]);
        end
    end
    
    if bDoZwicker
        for j = 1:length(h_Loudness)
            Saveas(h_Loudness(j),[paths.outputs 'zwicker-non-stat' num2str(j+1)]);
        end
        
        for j = 1:length(hZwicker_band)
            Saveas(hZwicker_band(j),[paths.outputs 'zwicker-per-band' num2str(j)]);
        end
        
    end
    
end

if bDiary
    diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end