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
% Last update on: 04/08/2014 % Update this date manually
% Last use on   : 04/08/2014 % Update this date manually
% 
% Original file name: VoD_comparisons.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 
    bDiary = 1;
end
Diary(mfilename,bDiary)

                        % tested on 01/07/2014
bDo_STFT        = 1;    % Analysis 1
bDo_load_audio  = bDo_STFT;    % Required to perform Analysis 1
bDoZwicker      = 0;    % Analysis replaced by bDoChalupper
bDoChalupper    = 1;
bMIR            = 0;    % Analysis discarded according to meeting on 23/07/2014

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
hMeas = [];
hModel = [];
hMeasr = [];
hModelr = [];
idx_i = 0; % to count figure handles

close all;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Signal alignment:
%       - You get the following plots (if info.bPlot == 1):
%           * 4 x Near and far field aligned signals (1 per acoustic mode)
if bDo_load_audio
    info.bSave = ~info.bSave; % to not to save again aligned files
    [ymeasured ymodelled misc] = VoD_read_aligned(bHPF,info);
    info.bSave = ~info.bSave;
else % only gets VoD params
    misc       = Get_VoD_params(1);
end
fs = misc.fs;
% Show_figures_one_by_one(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' type{idx} ' field']);
 
for mode = info.modes2check
    
    mode_idx = mode-1;
    
    if bDo_load_audio
        ymeas       = eval(['ymeasured.y' type{idx} num2str(mode) '_short;']);
        ymodel      = eval(['ymodelled.y' type{idx} num2str(mode) '_short;']);
    end
    
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
                tmp = Get_VoD_params(bGetT);
                t_coil_noalign = tmp.t_coil; % Correct this with the use of .t_coil
                t_coil = t_coil_noalign-repmat(tmp.ti_measured',1,length(t_coil_noalign));
                t_coil( find(t_coil < 0) ) = NaN;
                t_coil = Delete_NaN_columns(t_coil);
            end
            
            CurrT = misc.Tmodel(mode_idx);
            pos     = 0;% - misc.ti_model(mode_idx);
            
            Freq_res = misc.mf(mode_idx);
            
            subplot(2,1,1)
            hold on;
            plot(get(gca,'XLim'),[Freq_res Freq_res],'g--')
            
            if isfield(misc,'f0m1') % ACF values from Praat
                Exp1 = ['plot(misc.tf0' num2str(mode_idx) '/info.normalise_time_factor,misc.f0m' num2str(mode_idx) ',''r'')'];
                eval(Exp1)
            else
                % Predominant frequencies:
                fpmeas = Get_predominant_values(f,HdBmeas);
                plot(t/info.normalise_time_factor,fpmeas,'r') % centroid
            end

            subplot(2,1,2)
            hold on;
            plot(get(gca,'XLim'),[Freq_res Freq_res],'g--')
            
            if isfield(misc,'f0p1') % ACF values from Praat
                Exp1 = ['plot(misc.tf0' num2str(mode_idx) '/info.normalise_time_factor,misc.f0p' num2str(mode_idx) ',''r'')'];
                eval(Exp1)
            else
                % Predominant frequencies:
                fpmodel = Get_predominant_values(f,HdBmodel);
                plot(t/info.normalise_time_factor,fpmodel,'r') % centroid
            end

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
            h_STFT_Mesh(end+1) = gcf;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    % 3. Zwicker's model
    %       - Last used on: 03/07/2014
    if bDoZwicker
        
        VoD_comparisons_do_Zwicker;
         
    end
    
    % 4. Chalupper's model
    %       - Last used on: 01/08/2014
    
    if bDoChalupper
        ha1 = [];
        ha2 = [];    
        tmp_h = [];
        tmp_h = [];
        options = info;
        paths       = Get_TUe_subpaths('db_voice_of_dragon');
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filenames   = Get_filenames(paths.dir_calibrated_ms,'*.wav');
        
        options.nAnalyser = 12;
        [tmp_h tmp_ha]   = PsySoundCL([paths.dir_calibrated_ms filenames{mode_idx}],options); % 6 figures
        hMeas(10*idx_i+[1:6]) = tmp_h;
        ha1([1:6]) = tmp_ha;
        % 1. Average main loudness (Bark)       4. Main loudness (s)
        % 2. Average specific loudness (Bark)   5. Specific loudness (s)
        % 3. Loudness (s)                       6. Sharpness (s)
        
        options.nAnalyser = 15;
        [tmp_h tmp_ha]   = PsySoundCL([paths.dir_calibrated_ms filenames{mode_idx}],options); % 2 figures
        hMeas(10*idx_i+[7:8]) = tmp_h;
        ha1([7:8]) = tmp_ha;
        % 7. Specific roughness [Bark]
        % 8. Roughness [s]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filenames = Get_filenames(paths.dir_calibrated_ps,'*.wav');
        
        options.nAnalyser = 12;
        [tmp_h tmp_ha]  = PsySoundCL([paths.dir_calibrated_ps filenames{mode_idx}],options);
        hModel(10*idx_i+[1:6]) = tmp_h;
        ha2([1:6]) = tmp_ha;
        
        options.nAnalyser = 15;
        [tmp_h tmp_ha]  = PsySoundCL([paths.dir_calibrated_ps filenames{mode_idx}],options);
        hModel(10*idx_i+[7:8]) = tmp_h;
        ha2([7:8]) = tmp_ha;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i = 1:8
            
            YLimits = [get(ha1(i),'YLim'); get(ha2(i),'YLim')];
            YLimits = [min(min(YLimits)) max(max(YLimits))];
            linkaxes([ha1(i) ha2(i)],'y');
            set(ha1(i),'YLim',YLimits);
            
        end
        
        linkaxes([ha1([3:6,8]) ha2([3:6,8])],'x')
        xlim([0 3*misc.Tmodel(mode_idx)])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_i = idx_i+1;
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
    end
    
    if bDoZwicker
        for j = 1:length(h_Loudness)
            Saveas(h_Loudness(j),[paths.outputs 'zwicker-non-stat' num2str(j+1)]);
        end
        
        for j = 1:length(hZwicker_band)
            Saveas(hZwicker_band(j),[paths.outputs 'zwicker-per-band' num2str(j)]);
        end
        
    end
    
    if bDoChalupper
        nPlots = 8;
        idx = find(hMeas == 0);
        hMeas(idx) = [];
        idx = find(hModel == 0);
        hModel(idx) = [];
        ac_mode = ceil( hMeas /nPlots/2)+1;
        descrip = mod(hMeas,nPlots);
        
        for j = 1:length(hMeas)
            % figure(hMeas(j))
            option.bPrint = 0;
            txt = ['chalupper-ac' num2str(ac_mode(j)) '-descrip-' Num2str(descrip(j))];
            Saveas(hMeas(j),[paths.outputs  txt '-meas'],option);
            % figure(hModel(j));
            Saveas(hModel(j),[paths.outputs txt '-model'],option);
            
        end
        
    end
    
end

if bDiary
    diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end