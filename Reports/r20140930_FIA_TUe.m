function r20140930_FIA_TUe(options)
% function r20140930_FIA_TUe(options)
%
% 1. Description:
%       Generate figures and LaTeX tables for FIA 2014 paper (ID 0607). To
%       generate them, make sure you have changed manually bDoFluct, bDoAcMode2,
%       bDoAcMode5 to 1
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 0;
%       r20140930_FIA_TUe(options);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/09/2014
% Last update on: 26/09/2014 % Update this date manually
% Last use on   : 26/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all

bGenerate   = 0;
bDoFluct    = 0; % Figure 5
bDoAcMode5  = 1;
bDoAcMode2  = 0;
bDoAcMode2_filt = 0;

bDoDoppler  = 0;

bPlotOnly4paper = 1;

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
options = Ensure_field(options,'bSave',1);
options.bGenerateExcerpt = 0; % If 1, an excerpt is going to be generated

options.bDoAnalyser08 = 0; % SLM
options.bDoAnalyser10 = 0; %
options.bDoAnalyser11 = 1; % Fig 3, one-third OB analysis
options.bDoAnalyser12 = 1; % Loudness
options.bDoAnalyser15 = 1; % Roughness

bDiary = 0;
Diary(mfilename,bDiary);

%% Wav files used:
fi_ac5_meas     = [options.dest_folder 'meas-ac-mode-5.wav']; % measured
fi_ac5_model    = [options.dest_folder 'model-ac-mode-5.wav']; % predicted
fi_ac5_model_no_doppler = [options.dest_folder 'model-ac-mode-5-no-doppler.wav'];

fi_ac2_meas     = [options.dest_folder 'meas-ac-mode-2.wav']; % measured
fi_ac2_model    = [options.dest_folder 'model-ac-mode-2.wav']; % predicted
fi_ac2_model_no_doppler = [options.dest_folder 'model-ac-mode-2-no-doppler.wav'];

fi_ac2_meas_filt = [options.dest_folder 'meas-ac-mode-2-filt.wav']; % measured, filtered in Audacity

% Wavread(fi_ac2_meas);
% Wavread(fi_ac2_model);
% Wavread(fi_ac2_model_no_doppler);
% Wavread(fi_ac5_meas);
% Wavread(fi_ac5_model);
% Wavread(fi_ac5_model_no_doppler);

h = [];

tic

%% Plot options for Get_VoD_analysis:
options.LineWidth   = [1 2];
options.color       = {'b-','r--'};

%% Voice of the dragon params:

opts    = Get_VoD_params(0);
% Voice of the dragon
ceff    = 310;
L       = 0.7;
n       = 2:5;
fn      = n*ceff/(2*L);

%% new params

options.time2save   = 10;
options.N_periods2analyse = 3;

start_T_ac2 = 8;
ti2 = start_T_ac2 * opts.Tmodel(1);
tf2 = ti2+(1+options.N_periods2analyse)*opts.Tmodel(1);

start_T_ac5 = 3;
ti5 = start_T_ac5 * opts.Tmodel(4);
tf5 = ti5+(1+options.N_periods2analyse)*opts.Tmodel(4);

%% Table 1
tmp = [n' fn' opts.mf' ((fn-opts.mf)./fn*100)'];
numDecimals = 1;
tmp_col1_4 = Round(tmp,numDecimals);

tmp = [opts.Tmodel opts.Tmodel*4];
numDecimals = 4;
tmp_col5_6 = Round(tmp,numDecimals);
var2latex([tmp_col1_4 tmp_col5_6]);

%% Fluctuation strength, critical band levels LG (Figure 5)
if bDoFluct
    dBFS = 100; % Zwicker's calibration

    fi1 = fi_ac5_meas; % measured
    fi2 = fi_ac5_model; % predicted
    [x1 fs] = Wavread(fi1);
    x2 = Wavread(fi2);
    ttemp = ( 1:length(x1) )/fs;
    idx = find(ttemp > ti5-0.5 & ttemp < tf5+0.5 );
    x1 = x1(idx);
    x2 = x2(idx);
    
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    stPlot.Title1 = '(c) L_{Gmin}, ac. mode 5';
    stPlot.Title2 = '(d) L_{Gmax}, ac. mode 5';
    stPlot.Title3 = '(g)';
    stPlot.Title4 = '(h)';
    stPlot.Legend = {'meas','model'};
    stPlot.color = options.color;
    stPlot.LineWidth = options.LineWidth;
    [h(end+1:end+2) tmp] = Do_fluc_20140930(x1,x2,stPlot,fs);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Ac mode 5, fn is %.1f Hz, z(fn) is %.1f Bark\n',opts.mf(4),hz2bark(opts.mf(4)));
    idx = find(abs(tmp.diff_max)<=3 );
    disp('  - Frequencies where MAX are within 3 dB difference')
    disp(num2str(tmp.z(idx)));

    idx = find(abs(tmp.diff_min)<=3 );
    disp('  - Frequencies where MIN are within 3 dB difference')
    disp(num2str(tmp.z(idx)));

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %%%%
    fi1 = fi_ac2_meas; 
    fi2 = fi_ac2_model;
    [x1 fs] = Wavread(fi1);
    x2 = Wavread(fi2);
    ttemp = ( 1:length(x1) )/fs;
    idx = find(ttemp > ti2-0.5 & ttemp < tf2+0.5 );
    x1 = x1(idx);
    x2 = x2(idx);
    
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    stPlot.Title1 = '(a) L_{Gmin}, ac. mode 2';
    stPlot.Title2 = '(b) L_{Gmax}, ac. mode 2';
    stPlot.Title3 = '(e)';
    stPlot.Title4 = '(f)';
    [h(end+1:end+2) tmp] = Do_fluc_20140930(x1,x2,stPlot,fs);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Ac mode 2, fn is %.1f Hz, z(fn) is %.1f Bark\n',opts.mf(1),hz2bark(opts.mf(1)));
    idx = find(abs(tmp.diff_max)<=3 );
    disp('  - Frequencies where MAX are within 3 dB difference')
    disp(num2str(tmp.z(idx)));

    idx = find(abs(tmp.diff_min)<=3 );
    disp('  - Frequencies where MIN are within 3 dB difference')
    disp(num2str(tmp.z(idx)));

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    %%%%

    if bDoDoppler
        fi2 = fi_ac5_model_no_doppler;
        fi1 = fi_ac5_model;

        stPlot.YLim_fig1 = [0 62];
        stPlot.YLim_fig2 = [-5 27];
        stPlot.Legend = {'meas','model'};
        [h(end+1:end+2)] = Do_fluc_20140930(fi1,fi2,stPlot);
    end

    if options.bSave

        switch bPlotOnly4paper
            case 1
                g = Figure2paperfigure([h(1) h(2)],2);
                Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LG-ac5.eps'])
                g = Figure2paperfigure([h(3) h(4)],2);
                Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LG-ac2.eps'])
            case 0 
                Saveas(h(1),[options.dest_folder_fig 'psycho-fig-fluct-ac-5-min.eps'])
                Saveas(h(2),[options.dest_folder_fig 'psycho-fig-fluct-ac-5-max.eps'])
                Saveas(h(3),[options.dest_folder_fig 'psycho-fig-fluct-ac-2-min.eps'])
                Saveas(h(4),[options.dest_folder_fig 'psycho-fig-fluct-ac-2-max.eps'])
                Saveas(h(5),[options.dest_folder_fig 'psycho-fig-fluct-ac-5-max-doppler.eps'])
                Saveas(h(6),[options.dest_folder_fig 'psycho-fig-fluct-ac-5-max-doppler.eps'])
        end
    end
end

%% VoD analysis

options.calmethod   = 2; % 0 = AMT; 1 = dB(A)
options.bPlot       = 1;
options.modes2check = [2 5];
options.SPLrange = [10 90]; % for Figure 2
options.frange  = [50 5000]; % for Figure 2

options.ylim_loudness = [2.5 8];            % for Figure 4 (a,b), ylim
options.ylim_specific_loudness = [0 1.4];   % for Figure 4 (c,d), ylim
options.ylim_specific_roughness = [0 0.14]; % for Figure 6
options.ylim_bExtend = 1;

% Other 'options' fields:
%   - options.label1 - varies from plot to plot in PsySound
%   - options.label2 - varies from plot to plot in PsySound
%   - options.tanalysis - varies for each acoustic mode
%   - options.label  - EXPLAIN

% Options stft
overlap = 50;
nfft = 8192*2;
wlen = 2048; % 4096/10000 around 40 ms
nwtype = 4; % Hamming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of SPL as a function of loudness units
if bPlotOnly4paper == 0
    
    spl2lu(20:10:90);
    h = gcf;
    if options.bSave == 1
        Saveas(h, [options.dest_folder_fig 'fig-spl-LU']);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (Figure 2)
% Add window
if options.bSave == 1 & bGenerate == 0;
    options.bSave = ~options.bSave;
    out = VoD_write_aligned(options); % We do not want to generate again wav-files
    options.bSave = ~options.bSave;
else
    out = VoD_write_aligned(options);
end
 
set( out.ha(1),'XLim',[ti2 tf2]);
fprintf('ti = %.3f, tf = %.3f',ti2,tf2); % Use this to adapt Praat times
out.hFig(1)

ti5 = options.N_periods2analyse*out.Tmodel(4);
tf5 = ti5+(1+options.N_periods2analyse)*out.Tmodel(4);
set( out.ha(2),'XLim',[ti5 tf5]);
fprintf('ti = %.3f, tf = %.3f',ti5,tf5); % Use this to adapt Praat times

if options.bSave
    switch bPlotOnly4paper
        case 0
            for i = 1:length(out.hFig)
                options.format = 'epsc';
                Saveas(out.hFig(i)  ,[options.dest_folder_fig 'aligned-ac-mode-' num2str(options.modes2check(i))],options);
            end
        case 1
            g = Figure2paperfigure([out.hFig(1)],4);
            Saveas(g,[options.dest_folder_fig 'aligned-ac-mode-2'],options);
            g = Figure2paperfigure([out.hFig(2)],4);
            Saveas(g,[options.dest_folder_fig 'aligned-ac-mode-5'],options);
    end
end

%%

if bGenerate
    
    for i = 1:length(options.modes2check)
        
        lvl = 60;
        mode_idx = options.modes2check(i)-1;
        [x1 fs] = Wavread(out.Excerpt_m{mode_idx});
        x1 = Do_calibration_level(options.calmethod,x1,lvl);
        if options.bSave
            disp(['Overwrite with calibrated signal (across modes), adjusted to ' Num2str(lvl) 'dBSPL'])
            Wavwrite(x1,fs,out.Excerpt_m{mode_idx}); % Overwrite with calibrated signal
            [bSuccess bMessage] = movefile([out.Excerpt_m{mode_idx} '.wav'], options.dest_folder);
            if bSuccess == 0
                disp('File not moved, because: ')
                disp(bMessage)
            end
        end

        sil     = Gen_silence(2,fs);

        [x2 fs] = Wavread(out.Excerpt_p{mode_idx});
        x2 = Do_calibration_level(options.calmethod,x2,lvl); 
        if options.bSave
            Wavwrite(x2,fs,out.Excerpt_p{mode_idx});
            [bSuccess bMessage] = movefile([out.Excerpt_p{mode_idx} '.wav'], options.dest_folder);
            if bSuccess == 0
                disp('File not moved, because: ')
                disp(bMessage)
            end
        end

    end
    
    % File to be used in Praat
    
end

%% Figure 3, 4, 6, acoustic mode 2
if bDoAcMode5
    
    t1 = ti5;
    t2 = tf5;
    if bDoDoppler
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Acoustic mode 5, modelled with-without Doppler
        
        fi1 = fi_ac5_model_no_doppler; % predicted
        fi2 = fi_ac5_model; % predicted
        options.label1 = 'Model, no doppler,';
        options.label2 = 'Model';
        [x1 fs] = Wavread(fi1);
        [x2   ] = Wavread(fi2);

        if abs( rmsdb(x1) - rmsdb(x2) ) > 3
            yxfi1 = setdbspl(x1,rmsdb(x2)+100);
            Wavwrite(yxfi1,fs,fi2);
        end

        if bPlotOnly4paper == 0
            Get_LPC_frames(x1(fs:2*fs),fs);
            if options.bSave
                h = gcf;
                Saveas(h,[options.dest_folder_fig 'model-ac-meas-5-formants-no-doppler']);
            end
        end

        options.tanalysis = [t1 t2];
        options.label = 'ac-5-doppler';
        
        Get_VoD_analysis(fi1,fi2,options)
    end

    % Acoustic mode 5, measured and modelled
    fi1 = fi_ac5_meas; % measured
    fi2 = fi_ac5_model; % predicted
    options.label1 = ' meas';
    options.label2 = 'model';

    [x1 fs] = Wavread(fi1);
    [x2 fs] = Wavread(fi2);

    % STFT meas ac mode 5
    % rmsdbA(x1) % if -40 dB then is already calibrated to 60 dB(A)
    x1cal = x1;
    idx = find(To_dB(x1cal)<-100);
    x1cal(idx) = From_dB(-100);
    opts.scale_dB = 100; % 0 dBFS = 100 dB
    opts.txtTitle = ['Amplitude spectrogram, meas. hummer, ac. mode 5'];
    opts.XLabel = ['Time [s]'];
    figure;
    stft(x1cal, fs, nfft, wlen, overlap, nwtype,opts);
    ylim([0 2500])
    xlim([t1 t2])
    
    if bPlotOnly4paper == 0
        if options.bSave
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'meas-ac-5-stft']);
        end
    end

    if bPlotOnly4paper == 0
        Get_LPC_frames(x1(fs:2*fs),fs);
        if options.bSave
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'model-ac-meas-5-formants']);
        end
    end
    
    if bPlotOnly4paper == 0
        Get_LPC_frames(x2(fs:2*fs),fs);
        if options.bSave
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'model-ac-model-5-formants']);
        end
    end
    options.tanalysis = [t1 t2];
    options.label = 'ac-5';
    perc_ac_5 = Get_VoD_analysis(fi1,fi2,options);
end

%%
if bDoAcMode2
    t1 = ti2;
    t2 = tf2;
    if bDoDoppler
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Acoustic mode 2, modelled with-without Doppler
        fi1 = fi_ac2_model_no_doppler; 
        fi2 = fi_ac2_model; 
        options.label1 = 'Model, no doppler,';
        options.label2 = 'Model';
        [x1 fs] = Wavread(fi1);
        [x2   ] = Wavread(fi2);

        if abs( rmsdb(x1) - rmsdb(x2) ) > 3
            yxfi1 = setdbspl(x1,rmsdb(x2)+100);
            Wavwrite(yxfi1,fs,fi2);
        end

        if bPlotOnly4paper == 0
            Get_LPC_frames(x1(fs:2*fs),fs);
            if options.bSave
                h = gcf;
                Saveas(h,[options.dest_folder_fig 'model-ac-meas-2-formants-no-doppler']);
            end
        end

        options.tanalysis = [t1 t2];
        options.label = 'ac-2-doppler';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PsySound analysis:
        Get_VoD_analysis(fi1,fi2,options)% Acoustic mode 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
    if bDoAcMode2_filt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Acoustic mode 2
        fi1 = fi_ac2_meas; 
        fi2 = fi_ac2_meas_filt; 
        options.label1 = 'Meas';
        options.label2 = 'Meas, filt';
        [x1 fs] = Wavread(fi1);
        [x2   ] = Wavread(fi2);

        options.tanalysis = [t1 t2];
        options.label = 'ac-2-filtered';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PsySound analysis:
        out_ac_2_filtered = Get_VoD_analysis(fi1,fi2,options)% Acoustic mode 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acoustic mode 2

    fi1 = fi_ac2_meas; % measured
    fi2 = fi_ac2_model; % predicted
    options.label1 = ' meas';
    options.label2 = 'model';

    [x1 fs] = Wavread(fi1);
    [x2 fs] = Wavread(fi2);
    % [flpc1 out1] = Get_LPC(x1(fs:2*fs),fs);
    % [flpc2 out2] = Get_LPC(x2(fs:2*fs),fs);

    if bPlotOnly4paper == 0
        Get_LPC_frames(x1(fs:2*fs),fs);
        if options.bSave
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'model-ac-meas-2-formants']);
        end
    end
    
    % STFT with same params as ac mode 5
    x1cal = x1;
    idx = find(To_dB(x1cal)<-100);
    x1cal(idx) = From_dB(-100);
    opts.scale_dB = 100; % 0 dBFS = 100 dB
    opts.txtTitle = ['Amplitude spectrogram, meas. hummer, ac. mode 2'];
    opts.XLabel = ['Time [s]'];
    figure;
    stft(x1cal, fs, nfft, wlen, overlap, nwtype,opts);
    ylim([100 900])
    xlim([t1 t2])
    if options.bSave
        if bPlotOnly4paper == 0
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'meas-ac-2-stft']);
        end
    end

    if bPlotOnly4paper == 0
        Get_LPC_frames(x2(fs:2*fs),fs);
        if options.bSave
            h = gcf;
            Saveas(h,[options.dest_folder_fig 'model-ac-model-2-formants']);
        end
    end

    options.tanalysis = [t1 t2];
    options.label = 'ac-2';
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PsySound analysis:
    perc_ac_2 = Get_VoD_analysis(fi1,fi2,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Percentiles:

T = options.N_periods2analyse+1;

if bDoAcMode2
    idx2_L = find(perc_ac_2.Lt>=ti2 & perc_ac_2.Lt<=tf2);
    idx2_R = find(perc_ac_2.Rt>=ti2 & perc_ac_2.Rt<=tf2);
end

if bDoAcMode5
    idx5_L = find(perc_ac_5.Lt>=ti5 & perc_ac_5.Lt<=tf5);
    idx5_R = find(perc_ac_5.Rt>=ti5 & perc_ac_5.Rt<=tf5);
end

% TL5 = options.tanalysis(2)/(T-1) min(diff(perc_ac_2.Lt))
numDecimals = 2; % to round

col1 = [];
col2 = [];
col3 = [];
col4 = [];
col5 = [];
col6 = [];

if bDoAcMode2 == 1
    col1 = [    col1; ...
                mean(perc_ac_2.pL_in1.p5(start_T_ac2:start_T_ac2+T-1)); ...
                mean(perc_ac_2.pL_in2.p5(start_T_ac2:start_T_ac2+T-1))];
        
    col2 = [    col2;...
                std(perc_ac_2.pL_in1.p5(start_T_ac2:start_T_ac2+T-1)); ...
                std(perc_ac_2.pL_in2.p5(start_T_ac2:start_T_ac2+T-1))];
                
end

if bDoAcMode5 == 1
    col1 = [col1; ...
            mean(perc_ac_5.pL_in1.p5(start_T_ac5:start_T_ac5+T-1)); ...
            mean(perc_ac_5.pL_in2.p5(start_T_ac5:start_T_ac5+T-1)) ];
        
    col2 = [col2; ...
            std(perc_ac_5.pL_in1.p5(start_T_ac5:start_T_ac5+T-1)); ...
            std(perc_ac_5.pL_in2.p5(start_T_ac5:start_T_ac5+T-1)) ];
end
    
if bDoAcMode2_filt == 1
    col1 = [col1; ...
            mean(out_ac_2_filtered.pL_in2.p5(start_T_ac2:start_T_ac2+T-1))];
    col2 = [col2; ...
            std(out_ac_2_filtered.pL_in2.p5(start_T_ac2:start_T_ac2+T-1))];
end

if bDoAcMode2 == 1
    col3 = [    col3; ...
                mean(perc_ac_2.pL_in1.p50(start_T_ac2:start_T_ac2+T-1)); ...
                mean(perc_ac_2.pL_in2.p50(start_T_ac2:start_T_ac2+T-1))];

    col4 = [    col4; ...
                std(perc_ac_2.pL_in1.p50(start_T_ac2:start_T_ac2+T-1)); ...
                std(perc_ac_2.pL_in2.p50(start_T_ac2:start_T_ac2+T-1))];
end

if bDoAcMode5 == 1
    col3 = [    col3; ...
                mean(perc_ac_5.pL_in1.p50(start_T_ac5:start_T_ac5+T-1)); ...
                mean(perc_ac_5.pL_in2.p50(start_T_ac5:start_T_ac5+T-1)) ];

    col4 = [    col4; ...
                std(perc_ac_5.pL_in1.p50(start_T_ac5:start_T_ac5+T-1)); ...
                std(perc_ac_5.pL_in2.p50(start_T_ac5:start_T_ac5+T-1)) ];
end

if bDoAcMode2_filt == 1
    col3 = [col3; ...
            mean(out_ac_2_filtered.pL_in2.p50(start_T_ac2:start_T_ac2+T-1))];
    col4 = [col4; ...
            std(out_ac_2_filtered.pL_in2.p50(start_T_ac2:start_T_ac2+T-1))];
end

if bDoAcMode2 == 1
    col5 = [    col5; ...
                mean(perc_ac_2.pL_in1.p95(start_T_ac2:start_T_ac2+T-1)); ...
                mean(perc_ac_2.pL_in2.p95(start_T_ac2:start_T_ac2+T-1))];

    col6 = [    col6; ...
                std(perc_ac_2.pL_in1.p95(start_T_ac2:start_T_ac2+T-1)); ...
                std(perc_ac_2.pL_in2.p95(start_T_ac2:start_T_ac2+T-1))];
end

if bDoAcMode5 == 1
    col5 = [    col5; ...
                mean(perc_ac_5.pL_in1.p95(start_T_ac5:start_T_ac5+T-1)); ...
                mean(perc_ac_5.pL_in2.p95(start_T_ac5:start_T_ac5+T-1)) ];

    col6 = [    col6; ...
                std(perc_ac_5.pL_in1.p95(start_T_ac5:start_T_ac5+T-1)); ...
                std(perc_ac_5.pL_in2.p95(start_T_ac5:start_T_ac5+T-1)) ];
end


if bDoAcMode2_filt == 1
    col5 = [col5; ...
            mean(out_ac_2_filtered.pL_in2.p95(start_T_ac2:start_T_ac2+T-1))];
    col6 = [col6; ...
            std(out_ac_2_filtered.pL_in2.p95(start_T_ac2:start_T_ac2+T-1))];
end

% col7 = [    mean(perc_ac_2.pRmeas.p5(start_T_ac2:start_T_ac2+T-1)); ...
%             mean(perc_ac_2.pRmodel.p5(start_T_ac2:start_T_ac2+T-1)); ...
%             mean(perc_ac_5.pRmeas.p5(2:T+1)); ...
%             mean(perc_ac_5.pRmodel.p5(2:T+1)) ];
    
tmp = [];
for i = 1:6; 
    exp1 = ['tmp = [tmp col' num2str(i) '];'];
    eval(exp1);
end

if bDoAcMode2
    % Plots only to check table values
    figure;
    subplot(2,1,1)
    plot(   perc_ac_2.Lt(idx2_L), perc_ac_2.L_in1(idx2_L),...
            perc_ac_2.Lt(idx2_L), perc_ac_2.L_in2(idx2_L) )
    legend('meas','model')
    title('Loudness, as used for percentiles')
end

if bDoAcMode5
    subplot(2,1,2)
    plot(   perc_ac_5.Lt(idx5_L), perc_ac_2.L_in1(idx5_L),...
            perc_ac_5.Lt(idx5_L), perc_ac_2.L_in2(idx5_L) )
    legend('meas','model')
    title('Loudness, as used for percentiles')
end

%% Table 2
tmp = Round(tmp,numDecimals);
var2latex(tmp)

% Values from log:
%

% %  AC mode 5
% sum(DataLoud1): 
%   : 4.0502
% sum(DataLoud2): 
%   : 3.5548
% sum(DataLoud1-DataLoud2): 
%   : 0.49537

% %  AC mode 2:
% sum(DataLoud1): 
%   : 6.1561
% sum(DataLoud2): 
%   : 4.6768
% sum(DataLoud1-DataLoud2): 
%   : 1.4794

toc

%%
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
