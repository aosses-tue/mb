function r20141107_Perception_day_TUe(options)
% function r20141107_Perception_day_TUe(options)
%
% 1. Description:
%       Generate figures and LaTeX tables for FIA 2014 paper (ID 0607). To
%       generate them, make sure you have changed manually bDoFluct, bDoAcMode2,
%       bDoAcMode5 to 1
%
% 2. Stand-alone example:
%       options.bSave = 0;
%       r20141107_Perception_day_TUe(options);
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file : r20140930_FIA_TUe, r20141104_Acoustics_TUe 
% Created on    : 31/10/2014
% Last update on: 31/10/2014 % Update this date manually
% Last use on   : 31/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all

bDoFluct    = 0; % Not used for Dag van de perceptie
bDoFluct3   = 0;
bDoPsySound = 1;

options.bDoAnalyser01 = 1; % OK FFT
options.bDoAnalyser08 = 0; %    SLM
options.bDoAnalyser10 = 0; %
options.bDoAnalyser11 = 0; % Fig 3, one-third OB analysis
options.bDoAnalyser12 = 1; % Loudness
options.bDoAnalyser15 = 0; % NO Roughness

bDoAcMode5  = 1;
bDoAcMode2  = 1;
bDoAcMode2_filt = 1;

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
options = Ensure_field(options,'bSave',1);
options.bGenerateExcerpt = 0; % If 1, an excerpt is going to be generated

bDiary = 0;
Diary(mfilename,bDiary);

%% Wav files used:
fi_ac5_meas     = [options.dest_folder 'meas-ac-mode-5.wav']; % measured
fi_ac5_model    = [options.dest_folder 'model-ac-mode-5.wav']; % predicted

fi_ac2_meas     = [options.dest_folder 'meas-ac-mode-2.wav']; % measured
fi_ac2_model    = [options.dest_folder 'model-ac-mode-2.wav']; % predicted

fi_ac2_meas_filt = [options.dest_folder 'meas-ac-mode-2-filt.wav']; % measured, filtered in Audacity

h = [];

h2 = [];
h2name = [];

h5 = [];
h5name = [];

tic

%% Plot options for Get_VoD_analysis:
options.LineWidth   = [1 2 1];
options.color       = {'b-','r--','ko'};

%% Voice of the dragon params:

opts    = Get_VoD_params(0);
% ceff    = 310;
% L       = 0.7;
% n       = 2:5;
% fn      = n*ceff/(2*L);

%% new params

options.time2save   = 10;
options.N_periods2analyse = 3;

start_T_ac2 = 8;
ti2 = start_T_ac2 * opts.Tmodel(1);
tf2 = ti2+(options.N_periods2analyse)*opts.Tmodel(1);

start_T_ac5 = 3;
ti5 = start_T_ac5 * opts.Tmodel(4);
tf5 = ti5+(options.N_periods2analyse)*opts.Tmodel(4);

%% Fluctuation strength, critical band levels LG (Figure 5)
if bDoFluct
    
    tmp_opts.bDoFluct = bDoFluct;
    tmp_opts.bDoFluct3 = 0;
    
    tmp_opts.ti = ti5; % initial analysis time
    tmp_opts.tf = tf5;
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    
    stPlot.Title1 = '(c) L_{Gmin}, ac. mode 5';
    stPlot.Title2 = '(d) L_{Gmax}, ac. mode 5';
    stPlot.Title3 = '(g)';
    stPlot.Title4 = '(h)';
    
    stPlot.Legend = {'meas','model'};
    stPlot.color = options.color;
    stPlot.LineWidth = options.LineWidth;
    
    [h(end+1:end+2)] = Analyse_audio_classic_psychoacoustics(fi_ac5_meas,fi_ac5_model,[],tmp_opts,stPlot);
        
    %%%%
    
    tmp_opts.ti = ti2; % initial analysis time
    tmp_opts.tf = tf2;
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    stPlot.Title1 = '(a) L_{Gmin}, ac. mode 2';
    stPlot.Title2 = '(b) L_{Gmax}, ac. mode 2';
    stPlot.Title3 = '(e)';
    stPlot.Title4 = '(f)';
    [h(end+1:end+2)] = Analyse_audio_classic_psychoacoustics(fi_ac2_meas,fi_ac2_model,[],tmp_opts,stPlot);
    
    %%%%
    if bDoAcMode2_filt
        
        tmp_opts.ti = ti2; % initial analysis time
        tmp_opts.tf = tf2;
        stPlot.YLim_fig1 = [-2 62];
        stPlot.YLim_fig2 = [-5 27];
        stPlot.Title1 = '(a) L_{Gmin}, ac. mode 2, filt.';
        stPlot.Title2 = '(b) L_{Gmax}, ac. mode 2, filt.';
        stPlot.Title3 = '(e)';
        stPlot.Title4 = '(f)';
        stPlot.Legend = {'meas','model','meas filt'};
        % [h(end+1:end+2)] = Analyse_audio_classic_psychoacoustics(fi_ac2_meas_filt,fi_ac2_model,[],tmp_opts,stPlot);
        
        tmp_opts.bDoFluct = 0;
        tmp_opts.bDoFluct3 = 1;
        [h(end+1:end+2)] = Analyse_audio_classic_psychoacoustics(fi_ac2_meas,...
                                                                 fi_ac2_model,...
                                                                 fi_ac2_meas_filt, tmp_opts,stPlot);
        
    end
    
    if options.bSave

        g = Figure2paperfigure([h(1) h(2)],2);
        Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LG-ac5.eps'])
        g = Figure2paperfigure([h(3) h(4)],2);
        Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LG-ac2.eps'])
        if bDoAcMode2_filt
            g = Figure2paperfigure([h(5) h(6)],2);
            Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LG-ac2-filt.eps'])
        end        
        
    end
    
end

%%
if bDoFluct3
    
    tmp_opts.bDoFluct = 0;
    tmp_opts.bDoFluct3 = bDoFluct3;
    
    tmp_opts.ti = ti5; % initial analysis time
    tmp_opts.tf = tf5;
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    
    stPlot.Title1 = '(b) ac. mode 5';
    stPlot.Title2 = '';
    stPlot.Title3 = '';
    stPlot.Title4 = '';
    
    stPlot.Legend = {'meas','model'};
    stPlot.color = options.color;
    stPlot.LineWidth = options.LineWidth;
    
    [h(end+1)] = Analyse_audio_classic_psychoacoustics(fi_ac5_meas,fi_ac5_model,[],tmp_opts,stPlot);
        
    %%%%
          
    tmp_opts.ti = ti2; % initial analysis time
    tmp_opts.tf = tf2;
    stPlot.YLim_fig1 = [-2 62];
    stPlot.YLim_fig2 = [-5 27];
    stPlot.Title1 = '(a) ac. mode 2';
    stPlot.Title2 = '';
    stPlot.Title3 = '';
    stPlot.Title4 = '';
    stPlot.Legend = {'meas','model','meas filt'};

    tmp_opts.bDoFluct = 0;
    tmp_opts.bDoFluct3 = 1;
    [h(end+1)] = Analyse_audio_classic_psychoacoustics(fi_ac2_meas,...
                                                             fi_ac2_model,...
                                                             fi_ac2_meas_filt, tmp_opts,stPlot);
    
    if options.bSave

        g = Figure2paperfigure([h(end)],1,1);
        Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LGmax-2.eps'])        
        
        g = Figure2paperfigure([h(end-1)],1,1);
        Saveas(g,[options.dest_folder_fig 'psycho-fig-fluct-LGmax-5.eps'])        
        
    end
    
end

%% VoD analysis

options.calmethod   = 2; % 0 = AMT; 1 = dB(A)
options.bPlot       = 1;
options.modes2check = [2 5];
options.frange  = [50 5000]; % for Figure 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (Figure Waveforms)

stPlot = [];
stPlot.color        = options.color;
stPlot.LineWidth    = options.LineWidth;
stPlot.bPlotHorLine = 1;

acmode = 2;
stPlot.Title1 = sprintf('(a) meas, ac-mode %.0f',acmode);
stPlot.Title2 = sprintf('(c) modelled, ac-mode %.0f',acmode);
[h ha] = Get_VoD_waveform_2_5(fi_ac2_meas, fi_ac2_model, acmode, options,stPlot); % We do not want to generate again wav-files
 
T = opts.Tmodel(acmode-1);
set( ha(1),'XLim',[ti2 tf2]);
fprintf('ti = %.3f, tf = %.3f',ti2,tf2); 

acmode = 5;
stPlot.Title1 = sprintf('(b) meas, ac-mode %.0f',acmode);
stPlot.Title2 = sprintf('(d) modelled, ac-mode %.0f',acmode);
[h(end+1) ha(end+1:end+2)] = Get_VoD_waveform_2_5(fi_ac5_meas, fi_ac5_model, acmode, options,stPlot); % We do not want to generate again wav-files
T = opts.Tmodel(acmode-1);

set( ha(3),'XLim',[ti5 tf5]);
fprintf('ti = %.3f, tf = %.3f',ti5,tf5); % Use this to adapt Praat times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.bSave
    g1 = Figure2paperfigure(h(1),1);
    g2 = Figure2paperfigure(h(2),1);
    Saveas(g1,[options.dest_folder_fig 'Waveforms-ac-2'],options);
    Saveas(g2,[options.dest_folder_fig 'Waveforms-ac-5'],options);
    
    g  = Figure2paperfigure(h,2);
    
    Saveas(g,[options.dest_folder_fig 'Waveforms-ac-2-5'],options);
end

%% bGenerate eliminated from this script

%% Common params PsySound:

stPlot.ylim_loudness = [2.5 8];            % for Figure 4 (a,b), ylim
stPlot.ylim_specific_loudness = [0 1.5];   % for Figure 4 (c,d), ylim
stPlot.ylim_specific_roughness = [0 0.14]; % for Figure 6
stPlot.ylim_bExtend = 1;

%% Figure 3, 4, 6, acoustic mode 2
if bDoAcMode5
    
    acmode = 5;
    t1 = ti5;
    t2 = tf5;
    
    % Acoustic mode 5, measured and modelled
    fi1 = fi_ac5_meas; % measured
    fi2 = fi_ac5_model; % predicted
    stPlot.label1 = ' meas';
    stPlot.label2 = 'model';

    options.tanalysis = [t1 t2];
    options.label = 'ac-5';
    options.frange_FFT  = [950 1200];
    options.SPLrange    = [10 75];
    
    stPlot.Title = sprintf('ac-mode %.0f',acmode);          % FFT
    tmp_out = Get_PsySound_analysis(fi1,fi2,options,stPlot);% FFT
    
    h5 = [h5; tmp_out.h];
    h5name = [h5name; tmp_out.hname];
    
    options = rmfield(options,'frange_FFT');
    options = rmfield(options,'SPLrange');
end

%%
if bDoAcMode2
    acmode = 2;
    t1 = ti2;
    t2 = tf2;
    options.tanalysis = [t1 t2];
    options.ylim_roughness = [0 1.8];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acoustic mode 2

    fi1 = fi_ac2_meas; % measured
    fi2 = fi_ac2_model; % predicted
    fi3 = fi_ac2_meas_filt;
    
    stPlot.label1 = 'meas';
    stPlot.label2 = 'model';
    stPlot.label3 = 'meas,filt';

    options.label = 'ac-2';
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PsySound analysis:
    options.frange_FFT  = [300 600];    % FFT
    options.SPLrange    = [10 75];      % FFT
    
    stPlot.Title = sprintf('ac-mode %.0f',acmode);
    
    if bDoAcMode2_filt
        options.color{3} = 'k-';
        options.LineWidth(3) = 0.5;
        tmp_out = Get_PsySound_analysis_3_files(fi1,fi2,fi3,options,stPlot);
    else
        tmp_out = Get_PsySound_analysis(fi1,fi2,options,stPlot);
    end
    
    h2 = [h2; tmp_out.h];
    h2name = [h2name; tmp_out.hname];
    
    options = rmfield(options,'frange_FFT');
    options = rmfield(options,'SPLrange');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.bSave
        for i = 1:length(h2)
            
            figoption.format = 'epsc';
            
            g = Figure2paperfigure([h2(i)],1,1);
            Saveas(g,[options.dest_folder_fig h2name{i}],figoption);
            
            g = Figure2paperfigure([h5(i)],1,1);
            Saveas(g,[options.dest_folder_fig h5name{i}],figoption);
            
        end
    end
    
    % if bDoAcMode2_filt
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % Acoustic mode 2
    %     fi1 = fi_ac2_meas_filt; 
    %     fi2 = fi_ac2_model; 
    %     options.label1 = ' meas, filt';
    %     options.label2 = 'model';
    % 
    %     options.label = 'ac-2-filtered';
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % PsySound analysis:
    %     out_ac_2_filtered = Get_PsySound_analysis(fi1,fi2,options)% Acoustic mode 2
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end
end

toc

%%
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
