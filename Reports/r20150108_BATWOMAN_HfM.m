function r20150108_BATWOMAN_HfM(options)
% function r20150108_BATWOMAN_HfM(options)
%
% 1. Description:
%       Generate figures and LaTeX tables for FIA 2014 paper (ID 0607). To
%       generate them, make sure you have changed manually bDoFluct, bDoAcMode2,
%       bDoAcMode5 to 1
%
% 2. Stand-alone example:
%       options.bSave = 0;
%       r20150108_BATWOMAN_HfM(options);
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file : r20141104_Acoustics_TUe 
% Created on    : 07/01/2015
% Last update on: 07/01/2015 % Update this date manually
% Last use on   : 07/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

close all

bGenerate   = 0;
bDoFluct    = 0; % Figure 5
bDoAcMode5  = 1;
bDoAcMode2  = 1;
bDoAcMode2_filt = 1;
bCalculateRoughness = 0; % calculates roughness of existing test tones

bPlotOnly4paper = 1;

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio-HfM'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure-HfM' delim];
options = Ensure_field(options,'bSave',1);
options.bGenerateExcerpt = 0; % If 1, an excerpt is going to be generated

options.bDoAnalyser08 = 0; % SLM
options.bDoAnalyser10 = 0; %
options.bDoAnalyser11 = 0; % Fig 3, one-third OB analysis
options.bDoAnalyser12 = 1; % Loudness
options.bDoAnalyser15 = 1; % Roughness

bDiary = 0;
Diary(mfilename,bDiary);

%%
%%

close all
figure;
subplot(1,3,1)
plot([1 1], [0 3.2], 'k','LineWidth',2), hold on
grid on
deltax = 0.3; 
plot([1-deltax 1+deltax], [0.07   0.07] , 'g','LineWidth',.5), hold on
plot([1-deltax 1+deltax], [0.474  0.474], 'g','LineWidth',1), hold on
plot([1-deltax 1+deltax], [0.897  0.897], 'g','LineWidth',2), hold on
ylabel('Roughness [asper]')
set(gca,'XTickLabel','')
ylim([0 1])
Saveas(gcf, [options.dest_folder_fig 'R-1'])

plot([1-deltax 1+deltax], [0.01  0.01], 'r--','LineWidth',3), hold on % ac 2, model
plot([1-deltax 1+deltax], [0.21  0.21], 'b--','LineWidth',3), hold on % ac 2, meas
Saveas(gcf, [options.dest_folder_fig 'R-1-a2'])

plot([1-deltax 1+deltax], [0.15  0.15], 'r->','LineWidth',2), hold on % ac 5, model
plot([1-deltax 1+deltax], [0.13  0.13], 'b->','LineWidth',2), hold on % ac 5, meas
Saveas(gcf, [options.dest_folder_fig 'R-1-a2-5'])

%%

%% Wav files used:
fi_ac5_meas     = [options.dest_folder 'meas-ac-mode-5.wav']; % measured
fi_ac5_model    = [options.dest_folder 'model-ac-mode-5.wav']; % predicted

fi_ac2_meas     = [options.dest_folder 'meas-ac-mode-2.wav']; % measured
fi_ac2_model    = [options.dest_folder 'model-ac-mode-2.wav']; % predicted

fi_ac2_meas_filt = [options.dest_folder 'meas-ac-mode-2-filt.wav']; % measured, filtered in Audacity

h = [];

tic

%% Plot options for Get_MI_analysis:
options.LineWidth   = [1 2];
options.color       = {'b-','r--'};

%% Voice of the dragon params:

opts    = Get_VoD_params(0);
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

if bCalculateRoughness
    local_dir = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20150108-10-BATWOMAN-HfM\Audio\';
    test_files = {  'rough_ref.wav',...
                    'test_fc_1500_AM_m_018_fmod_070Hz.wav', ...
                    'test_fc_1500_AM_m_053_fmod_070Hz.wav', ...
                    'test_fc_1500_AM_m_090_fmod_070Hz.wav'};
    tmp = [];
	tmp.nAnalyser = 15;
    tmp.CalMethod = 2; % Zwicker
    
    tmp.ti = 0.2;
    tmp.tf = 0.8;
	for i = 1:3
        tmp_out = PsySoundCL([local_dir test_files{i}], tmp);
        disp('')
        
        idx = find( tmp_out.t>=tmp.ti & tmp_out.t<=tmp.tf );
        tmp_res = mean( tmp_out.DataRough(idx));
        fprintf('test file %s = %.3f [asper]\n\n',test_files{i},tmp_res);
        
    end
end

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

    if bDoAcMode2_filt
        %%%%
        fi1 = fi_ac2_meas_filt; 
        fi2 = fi_ac2_model;
        [x1 fs] = Wavread(fi1);
        x2 = Wavread(fi2);
        ttemp = ( 1:length(x1) )/fs;
        idx = find(ttemp > ti2-0.5 & ttemp < tf2+0.5 );
        x1 = x1(idx);
        x2 = x2(idx);

        stPlot.YLim_fig1 = [-2 62];
        stPlot.YLim_fig2 = [-5 27];
        stPlot.Title1 = '(a) L_{Gmin}, ac. mode 2, filt.';
        stPlot.Title2 = '(b) L_{Gmax}, ac. mode 2, filt.';
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
    end
    if options.bSave

        switch bPlotOnly4paper
            case 1
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
    
end

%% Figure 3, 4, 6, acoustic mode 2
if bDoAcMode5
    
    t1 = ti5;
    t2 = tf5;
    
    % Acoustic mode 5, measured and modelled
    fi1 = fi_ac5_meas; % measured
    fi2 = fi_ac5_model; % predicted
    options.label1 = ' meas';
    options.label2 = 'model';

    options.tanalysis = [t1 t2];
    options.label = 'ac-5';
    perc_ac_5 = Get_MI_analysis(fi1,fi2,options);
    
    idx = find(perc_ac_5.Rt>=options.tanalysis(1) & perc_ac_5.Rt<=options.tanalysis(2));
    ro1 = mean( perc_ac_5.R_in1(idx) );
    ro2 = mean( perc_ac_5.R_in2(idx) );
    fprintf('Ravg 1 meas, ac mode 5 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro1);
    fprintf('Ravg 2 mod , ac mode 5 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro2);
    
end

%%
if bDoAcMode2
    t1 = ti2;
    t2 = tf2;
    options.tanalysis = [t1 t2];
    options.ylim_roughness = [0 1.8];
    
    if bDoAcMode2_filt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Acoustic mode 2
        fi1 = fi_ac2_meas_filt; 
        fi2 = fi_ac2_model; 
        options.label1 = ' meas, filt';
        options.label2 = 'model';

        options.label = 'ac-2-filtered';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PsySound analysis:
        options.ylim_bExtend = 1;
        options.ylim_bDrawLine = 0;
        out_ac_2_filtered = Get_MI_analysis(fi1,fi2,options);% Acoustic mode 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %idx = find(perc_ac_5.Rt        >=options.tanalysis(1) & perc_ac_5.Rt<=options.tanalysis(2));
        idx = find(out_ac_2_filtered.Rt>=options.tanalysis(1) & out_ac_2_filtered.Rt<=options.tanalysis(2));
        ro1 = mean( out_ac_2_filtered.R_in1(idx) );
        ro2 = mean( out_ac_2_filtered.R_in2(idx) );
        fprintf('Ravg 1 meas filt, ac mode 2 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro1);
        fprintf('Ravg 2 mod      , ac mode 2 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro2);
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acoustic mode 2

    fi1 = fi_ac2_meas; % measured
    fi2 = fi_ac2_model; % predicted
    options.label1 = ' meas';
    options.label2 = 'model';

    options.label = 'ac-2';
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PsySound analysis:
    perc_ac_2 = Get_MI_analysis(fi1,fi2,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    idx = find(perc_ac_2.Rt>=options.tanalysis(1) & perc_ac_2.Rt<=options.tanalysis(2));
    ro1 = mean( perc_ac_2.R_in1(idx) );
    ro2 = mean( perc_ac_2.R_in2(idx) );
    fprintf('Ravg 1 meas, ac mode 2 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro1);
    fprintf('Ravg 2 mod , ac mode 2 for t(%.2f,%.2f) = %.3f [asper] \n',t1,t2,ro2);
    
end

toc

%%
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
