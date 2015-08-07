function $$functionname$$(f1,f2)
% function $$functionname$$(f1,f2)
%
% 1. Description:
%       Reproduce processing done using PsySoundControl GUI.
% 
% 2. Stand-alone example: 
%       % This is supposed to be a stand-alone script so just run:
%       $$functionname$$;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : $$dd$$/$$mm$$/$$yyyy$$ (automatic generation)
% Last update on: $$dd$$/$$mm$$/$$yyyy$$ 
% Last use on   : $$dd$$/$$mm$$/$$yyyy$$ 
% Template      : template_PsySoundCL.m (revised on: 12/03/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% Copy lines 100 to 460

close all
if nargin < 1
    f1 = '$$f1$$';
end
if nargin < 2
    f2 = '$$f2$$';
end
dir_output = Get_TUe_paths('outputs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L102 on 12/03/2015
bUsePsySound = $$bUsePsySound$$;

options.bUsePsySound = bUsePsySound;

bSave       = $$bSave$$;
nAnalyser   = $$nAnalyser$$;

nSkipStart  = $$nSkipStart$$;
nSkipEnd    = $$nSkipEnd$$;

options.bLogScale = 0; % get(handles.cbLogAxis,'value'); % to be used in PsySound_Figures

fs          = $$fs$$;

sample_inf = $$sample_inf$$;
sample_sup = $$sample_sup$$;

HopSize = [$$HopSize$$];
CParams.HopSize = HopSize;

toffset = $$toffset$$; % time offset of audio 2 in relation to audio 1

bGenerateExcerpt = $$bGenerateExcerpt$$;
tanalysis = [$$tanalysis$$];

options.bGenerateExcerpt = bGenerateExcerpt;
options.tanalysis = tanalysis;

if bGenerateExcerpt
    options.tanalysis = [sample_inf-1 sample_sup-1]/fs;
    
  
    
end

eval( sprintf('options.bDoAnalyser%s=1;',Num2str($$nAnalyser$$,2)) ); % activates selected processor

filename1 = f1; %get(handles.txtFile1,'string');
filename2 = f2; %get(handles.txtFile2,'string');

options.nAnalyser   = nAnalyser;
options.bSave       = bSave;

if options.bSave == 1
    options.dest_folder_fig = dir_output;
end
    
%% Plot options:
options.label1 = '$$label1$$';
options.label2 = '$$label2$$';

options          = ef(options,'label','');
options.SPLrange = [$$SPLrange$$];
options.frange   = [$$frange$$];
options.zrange   = [$$zrange$$];
options.trange   = [$$trange$$];
options          = ef(options,'ylim_bExtend',0);
options          = ef(options,'ylim_bDrawLine',0);
options.bLoudnessContrained = $$bLoudnessContrained$$; % only validated for Analyser 12
options.zlim4assessment = options.zrange; 

if bUsePsySound
    
    options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
    
    callevel = $$callevel$$; % rms 90 dB SPL = 0 dBFS 

    options     = ef(options,'bDoAnalyser08',0);
    options     = ef(options,'bDoAnalyser10',0);
    options     = ef(options,'bDoAnalyser11',0);
    options     = ef(options,'bDoAnalyser12',0);
    options     = ef(options,'bDoAnalyser15',0);

    tmp_h = [];
    tmp_h = [];

    options = Ensure_field(options,'calfile',[Get_TUe_paths('db_calfiles') 'track_03.wav']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    options.bCosineRamp = 0; % Cos ramp not yet applied for Loudness calculations
    options.bCosineRampOnset = 0; %ms

    options.bPlot = 0;

    G1 = $$G1$$;
    G2 = $$G2$$;
    options.callevel = callevel + G1;
    options.nSkipStart = nSkipStart;
    [out_1 tmp_h tmp_ha]   = PsySoundCL(filename1,options,CParams);
    
    options.callevel = callevel + G2;
    options.nSkipStart = nSkipStart;
    [out_2 tmp_h tmp_ha]   = PsySoundCL(filename2,options,CParams);
    
else
    
    callevel = str2num( get(handles.txtCalLevel,'string') ); % rms 90 dB SPL = 0 dBFS 
    warning('callevel not used at all for m-file scripts...');
    
    [insig1 fs1] = Wavread(filename1);
    [insig2 fs2] = Wavread(filename2);
    
    if options.bGenerateExcerpt
        Nextra = length(insig1)-(sample_sup - sample_inf)-1;
        if Nextra >= 8192
            Nextra = 8192;
        end
            
        try
            insig1tmp = insig1( sample_inf:sample_sup + Nextra ); % one additional frame
            insig2tmp = insig2( sample_inf + toffset:sample_sup + toffset + Nextra ); % one additional frame
            insig1 = insig1tmp;
            insig2 = insig2tmp;
        catch
            insig1 = insig1( sample_inf:sample_sup ); % one additional frame
            insig2 = insig2( sample_inf + toffset:sample_sup + toffset ); % one additional frame
            warning('using catch...')
        end
        set(handles.txtExcerpt,'visible','on');
         
        %% Extracts excerpt:
        fname1_excerpt = [Delete_extension(filename1,'wav') '-e.wav']; 
        fname2_excerpt = [Delete_extension(filename2,'wav') '-e.wav']; 
 
        Wavwrite(insig1,fs1,fname1_excerpt);
        Wavwrite(insig2,fs2,fname2_excerpt);
        
        filename1 = fname1_excerpt;
        filename2 = fname2_excerpt;
        
    else
        set(handles.txtExcerpt,'visible','off');
    end
    
    lvl_m_30dBFS = str2num( get(handles.txtCalLevel,'string') );
    calvalue = lvl_m_30dBFS-60; % values calibrated to 90 dB RMS = 0 dBFS (Fastl's standard)
    insig1 = From_dB(calvalue+handles.audio.G1) * insig1;    
    insig2 = From_dB(calvalue+handles.audio.G2) * insig2;
    
    switch options.nAnalyser
        
        case 12 % Loudness
            
            % Only loudness fluctuation:
            dBFS = lvl_m_30dBFS + 30 - calvalue;
            
            tinsig = ( 0:length(insig1)-1 )/fs1;
            idx = find( tinsig>=options.trange(1) & tinsig<=options.trange(2) );
            
            if length(idx) ~= length(tinsig)
                insig1 = insig1(idx);
                insig2 = insig2(idx);
                warning('Truncating insig, check if it is working properly when using a non-zero off-set');
            end
                
            [xx out_1] = LoudnessFluctuation_offline(insig1,[],fs,dBFS);
            [xx out_2] = LoudnessFluctuation_offline(insig2,[],fs,dBFS);
              
        case 15 % Roughness
            
            N = 8192; % default frame length
            opts.nSkipStart = nSkipStart;
            [xx out_1] = Roughness_offline(insig1,fs1,N,opts,CParams,0);
            [xx out_2] = Roughness_offline(insig2,fs2,N,opts,CParams,0);
            Ndel = length(out_1.t);
            out_1.t(Ndel)       = []; % we delete last frame
            out_1.Data1(Ndel)   = [];
            out_1.Data2(Ndel,:) = [];
            out_2.t(Ndel)       = [];
            out_2.Data1(Ndel)   = [];
            out_2.Data2(Ndel,:) = [];
            
        case 20 % Fluctuation strength, see also r20141126_fluctuation
            
            N = 44100*4; % 8192*4;
            opts.nSkipStart = nSkipStart;
            warning('Fluctuation strength: temporal value...')
            [xx out_1] = FluctuationStrength_offline_debug(insig1(1:N),fs1,N,0);
            [xx out_2] = FluctuationStrength_offline_debug(insig2(1:N),fs2,N,0);
            
        case 21 % Calculation made at plot section
            warning('color of the series are not automated')
    end
    
end

param   = [];
h       = []; % handles figures
ha      = [];

bPlotParam1 = $$bPlotParam1$$; 
bPlotParam2 = $$bPlotParam2$$; 
bPlotParam3 = $$bPlotParam3$$; 
bPlotParam4 = $$bPlotParam4$$; 
bPlotParam5 = $$bPlotParam5$$; 
bPlotParam6 = $$bPlotParam6$$; 
bPlotParam7 = $$bPlotParam7$$; 

labelParam1 = '$$labelParam1$$';
labelParam2 = '$$labelParam2$$';
labelParam3 = '$$labelParam3$$';
labelParam4 = '$$labelParam4$$';
labelParam5 = '$$labelParam5$$';
labelParam6 = '$$labelParam6$$';
labelParam7 = '$$labelParam7$$';

%% Plots
if nAnalyser ~= 21
    
    bPercentiles = $$bPercentiles$$; % only for loudness
    
    if bPercentiles & options.nAnalyser == 12
        param{end+1} = 'loudness-percentiles';
        [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
        param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    end
    
    if bPlotParam1
        % Loudness, Roughness
        param{end+1}        = labelParam1;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
        param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

        if isfield(options,'ylim') % ylim_loudness, ylim_roughness
            ylim(options.ylim);
        end
        if isfield(options,'ylim_bExtend')

            if options.ylim_bExtend == 1
                y_old = get(gca,'YLim');
                x_old = get(gca,'XLim');
                ylim_extend(gca,1.25);

                if options.ylim_bDrawLine == 1
                    plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
                end
            end
        end

    end

    if bPlotParam2
        % Specific loudness, roughness
        param{end+1}        = labelParam2;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end

    if bPlotParam3
        param{end+1}        = labelParam3;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end

    if bPlotParam4
        param{end+1}        = labelParam4;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end

    if bPlotParam5
        param{end+1}        = labelParam5;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end

    if bPlotParam6
        param{end+1}        = labelParam6;
        [h(end+1) ha(end+1) stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end
    
    if bPlotParam7
        % 12 - Loudness fluctuation
        if strcmp(labelParam7,'loudness-fluctuation')
            param{end+1} = [labelParam7 '-max'];
            [h(end+1) xx  stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
            param{end+1} = [labelParam7 '-min'];
            [h(end+1) xx  stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
        else
            param{end+1}        = labelParam7;
            [h(end+1) xx  stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
        end
    end
    
else
    if nAnalyser == 21
        
        options.type = 4; % 40-ms frame length
        % Praat Analyser:
        [h(end+1:end+2)] = Get_waveforms_and_F0_praat(filename1,filename2,options);
        param{1} = 'time-series';
        
        param{2} = 'fundamental-frequency';
        
    end
end

param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assignin('base', 'h', h);
assignin('base', 'ha', ha);

if options.bSave
    
    disp('Figures are going to be stored... Press ctr+C to cancel')
    
    try
        paths.outputs = options.dst_folder_fig;
    catch
        paths.outputs   = Get_TUe_paths('outputs');
    end
    
    for i = 1:length(h)
        % options.format = 'emf';
        % Saveas(h(i),[options.dest_folder_fig 'psycho-fig-' param{i} options.label],options);
        options.format = 'fig';
        Saveas(h(i),[options.dest_folder_fig 'fig-' param{i} options.label],options);
        options.format = 'epsc';
        Saveas(h(i),[options.dest_folder_fig 'fig-' param{i} options.label],options);
        
        % Generating txt/log-file, once per analyser:
        if i == 1
            fndiary = [options.dest_folder_fig 'fig-log-analyser-' num2str(nAnalyser) '.txt'];
            diary(fndiary)

            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf('Date/time of processing: %s\n', Get_date_ddmmyyyy(1));
            fprintf('Output directory: %s\n'    ,options.dest_folder_fig);
            fprintf('Level ref. tone: %.1f dB\n',callevel);
            fprintf('File name 1: %s (gain = %.2f dB)\n',filename1,G1);
            fprintf('File name 2: %s (gain = %.2f dB)\n',filename2,G2);
            fprintf('Initial/final sample: %.0f, %.0f\n',sample_inf,sample_sup);
            fprintf('trange = (%.3f, %.3f) [s]\n',options.trange(1),options.trange(2));
            fprintf('label 1: %s',options.label1);
            fprintf('label 2: %s',options.label2);
            
            if bUsePsySound 
                fprintf('Processed using PsySound \n');
            else
                fprintf('Processed using m-files (not PsySound) \n');
            end
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
                
            diary off
        end
        % end Generating txt file
    end
    
else
    
    disp('Figures are NOT going to be stored... Set options.bSave to 1 and re-run the scripts in case you want to save the figures')
    pause(1)
    
end

if ~bUsePsySound
    
    disp('...deleting temporal audio file');
    delete( fname1_excerpt );
    delete( fname2_excerpt );
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
