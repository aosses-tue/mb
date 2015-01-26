function varargout = PsySoundControl(varargin)
% PSYSOUNDCONTROL MATLAB code for PsySoundControl.fig
%      PSYSOUNDCONTROL, by itself, creates a new PSYSOUNDCONTROL or raises the existing
%      singleton*.
%
%      H = PSYSOUNDCONTROL returns the handle to a new PSYSOUNDCONTROL or the handle to
%      the existing singleton*.
%
%      PSYSOUNDCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSYSOUNDCONTROL.M with the given input arguments.
%
%      PSYSOUNDCONTROL('Property','Value',...) creates a new PSYSOUNDCONTROL or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PsySoundControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PsySoundControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%   
%   Limitations (written by AO):
%       - Roughness calculations using PsySound include first and last frame
%       in the average value calculations. Use Roughness_offline.m
% 
%       Line    Stage                               Last updated on
%       37      Initialisation                      18/01/2015
%       94      Calculation - calculate_Callback    19/01/2015
%       267     reset_Callback                      18/01/2015
%       309     Initialisation GUI                  18/01/2015
%       313     unitgroup_SelectionChangeFcn        21/01/2015
%       479     popAnalyser_Callback                25/01/2015
%       783     Load data                           21/01/2015
%       
% TO DO:
%       2. SLM: problem at @Analyser/process, line 381. Object subclass is not readable 'AZ'
%       3. Put save figures in other Button
%       4. xlim axis, ylim axis
% 
% Edit the above text to modify the response to help PsySoundControl
% Created on        : 16/01/2015
% Last modified on  : 21/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PsySoundControl_OpeningFcn, ...
                   'gui_OutputFcn',  @PsySoundControl_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before PsySoundControl is made visible.
function PsySoundControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PsySoundControl (see VARARGIN)

% Choose default command line output for PsySoundControl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes PsySoundControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PsySoundControl_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Calculations: 
% Executes on button press in calculate.

function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Needed: 
%       Num2str
%       ef

bUsePsySound = get(handles.rbPsySound,'value');
bSave       = get(handles.bSave,'value');
nAnalyser   = get(handles.popAnalyser,'value');
dir_output  = get(handles.txtOutputDir,'string');

nSkipStart  = str2num( get(handles.txtAnalysisStart,'string') );
nSkipEnd    = str2num( get(handles.txtAnalysisEnd  ,'string') );

fs          = handles.audio.fs;
tanalysis_inf = str2num(get(handles.txtti,'string')); % in samples
tanalysis_sup = str2num(get(handles.txttf,'string')); % in samples

toffset = handles.audio.toffset; % time offset of audio 2 in relation to audio 1

options.bGenerateExcerpt = handles.audio.bGenerateExcerpt;

if options.bGenerateExcerpt
    options.tanalysis = [tanalysis_inf tanalysis_sup]/fs;
    set( handles.txtExcerpt,'visible','on');
 else
    set( handles.txtExcerpt,'visible','off');
end

eval( sprintf('options.bDoAnalyser%s=1;',Num2str(nAnalyser,2)) ); % activates selected processor

filename1 = get(handles.txtFile1,'string');
filename2 = get(handles.txtFile2,'string');

options.nAnalyser   = nAnalyser;
options.bSave       = bSave;

if options.bSave == 1
    options.dest_folder_fig = dir_output;
end
    
%% Plot options:
options.label1 = get(handles.txtLabel1,'string');
options.label2 = get(handles.txtLabel2,'string');

options     = Ensure_field(options,'bPlot',1);
options     = Ensure_field(options,'label','');
options     = Ensure_field(options,'SPLrange',[10 70]);
options.frange = [str2num(get(handles.txtFreqmin,'string')) str2num(get(handles.txtFreqmax,'string'))];
options     = Ensure_field(options,'ylim_bExtend',0);
options     = Ensure_field(options,'ylim_bDrawLine',0);
    
if bUsePsySound
    
    options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
    
    callevel = str2num( get(handles.txtCalLevel,'string') ); % rms 90 dB SPL = 0 dBFS 

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

    options.callevel = callevel + handles.audio.G1;
    options.nSkipStart = nSkipStart;
    [out_1 tmp_h tmp_ha]   = PsySoundCL(filename1,options);
    
    options.callevel = callevel + handles.audio.G2;
    options.nSkipStart = nSkipStart;
    [out_2 tmp_h tmp_ha]   = PsySoundCL(filename2,options);
    
else
    
    [insig1 fs1] = Wavread(filename1);
    [insig2 fs2] = Wavread(filename2);
    
    if options.bGenerateExcerpt
        Nextra = length(insig1)-(tanalysis_sup - tanalysis_inf)-1;
        if Nextra >= 8192
            Nextra = 8192;
        end
            
        insig1 = insig1( tanalysis_inf:tanalysis_sup + Nextra ); % one additional frame
        insig2 = insig2( tanalysis_inf + toffset:tanalysis_sup + toffset + Nextra ); % one additional frame
        set(handles.txtExcerpt,'visible','on');
    else
        set(handles.txtExcerpt,'visible','off');
    end
    
    calvalue = str2num( get(handles.txtCalLevel,'string') )-60; % values calibrated to 90 dB RMS = 0 dBFS (Fastl's standard)
    insig1 = From_dB(calvalue+handles.audio.G1) * insig1;    
    insig2 = From_dB(calvalue+handles.audio.G2) * insig2;
    
    switch options.nAnalyser
        
        case 12 % Loudness
            
            % Only loudness fluctuation:
            [xx out_1] = LoudnessFluctuation_offline(insig1,[],fs);
            [xx out_2] = LoudnessFluctuation_offline(insig2,[],fs);
              
        case 15 % Roughness
            
            N = 8192; % default frame length
            opts.nSkipStart = nSkipStart;
            [xx out_1] = Roughness_offline(insig1,fs1,N,opts,0);
            [xx out_2] = Roughness_offline(insig2,fs2,N,opts,0);
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

for i = 1:7
    exp1 = sprintf('bPlotParam%.0f = get(handles.chParam%.0f,''value''); labelParam%.0f = get(handles.chParam%.0f,''string'');',i,i,i,i);
    eval( exp1 );
end

%% Plots
if nAnalyser ~= 21
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
        param{end+1}        = labelParam7;
        [h(end+1) xx        stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    end
    
else
    if nAnalyser == 21
        
        %title( 'f0 extraction (Praat Analyser)' );
        [h ha] = Get_waveforms_and_F0_praat(filename1,filename2,options);
        param{1} = 'fundamental-frequency';
        
    end
end

param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assignin('base', 'h', h);
assignin('base', 'ha', ha);

if options.bSave
    
    disp('Figures are going to be stored... Press ctr+C to cancel')
    pause(2)
    
    try
        paths.outputs = options.dst_folder_fig;
    catch
        paths.outputs   = Get_TUe_paths('outputs');
    end
    
    for i = 1:length(h)
        % options.format = 'emf';
        % Saveas(h(i),[options.dest_folder_fig 'psycho-fig-' param{i} options.label],options);
        options.format = 'epsc';
        Saveas(h(i),[options.dest_folder_fig 'psycho-fig-' param{i} options.label],options);
    end
    
else
    
    disp('Figures are NOT going to be stored... Set options.bSave to 1 and re-run the scripts in case you want to save the figures')
    pause(2)
    
end

%% reset - Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

%% radio group change 
%       - Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nAnalyser = get(handles.popAnalyser,'value');

if (hObject == handles.rbPsySound)
    set(handles.txtScript, 'String', '');
else
    switch nAnalyser
        case 15
            set(handles.txtScript, 'String', 'Roughness_offline.m');
        case 20 
            set(handles.txtScript, 'String', 'FluctuationStrength_offline_debug.m');
        otherwise
            set(handles.txtScript, 'String', '');
    end
end

%% Initialisation GUI
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'audio') && ~isreset
    return;
end

handles.audio.fs = 0;
handles.audio.bGenerateExcerpt = 0;

set(handles.unitgroup, 'SelectedObject', handles.rbPsySound);
set(handles.txtExcerpt,'visible','off');

try
    set(handles.txtOutputDir,'string',Get_TUe_paths('outputs'));
catch
    warning('Type your output dir in the GUI');
end

% Update handles structure
guidata(handles.figure1, handles);

%% Checkboxes:
% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles)
% hObject    handle to bSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSave
bSave = get(handles.bSave,'value');

if bSave
    set(handles.txtOutputDir,'enable','on');
else
    set(handles.txtOutputDir,'enable','off');
end

% --- Executes on button press in chFile1.
function chFile1_Callback(hObject, eventdata, handles)
% hObject    handle to chFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bFile1 = get(handles.chFile1,'value');

if bFile1
    set(handles.txtFile1,'enable','on');
else
    set(handles.txtFile1,'enable','off');
end

% --- Executes on button press in chFile2.
function chFile2_Callback(hObject, eventdata, handles)
% hObject    handle to chFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chFile2
bFile2 = get(handles.chFile2,'value');

if bFile2
    set(handles.txtFile2,'enable','on');
else
    set(handles.txtFile2,'enable','off');
end

function chParam1_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chParam2_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chParam3_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chParam4_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chParam5_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chParam6_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% --- Executes on selection change in popAnalyser.
function popAnalyser_Callback(hObject, eventdata, handles)
% hObject    handle to popAnalyser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nAnalyser = get(handles.popAnalyser,'value');

switch nAnalyser

    case 1
        
        set(handles.txtAnalysisStart,'Enable','off');
        
        set(handles.rbPsySound,'value',1);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',0);
        set(handles.rbScripts,'enable','off');
        
        set(handles.chParam1,'string','spectrogram'); % 381 x 1024 x 381
        set(handles.chParam1,'enable','on');
        set(handles.chParam1,'value',0);
        
        set(handles.chParam2,'string','average-power-spectrum'); % 1 x 1024
        set(handles.chParam2,'enable','on');
        % set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','spectral-centroid'); % 
        set(handles.chParam3,'enable','off');
        set(handles.chParam3,'value',0);
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
    
    case 10
        
        set(handles.txtAnalysisStart,'Enable','off');
        
        set(handles.rbPsySound,'value',1);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',0);
        set(handles.rbScripts,'enable','off');
        
        set(handles.chParam1,'string','one-third-octave-band-spectrogram'); % 440 x 28 x 440
        set(handles.chParam1,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam2,'string','one-third-octave-band-spectrum'); % 1 x 28
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','specific-loudness'); % specific-loudness 1 x 240
        set(handles.chParam3,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
    
    case 11
        
        set(handles.txtAnalysisStart,'Enable','off');
        
        set(handles.rbPsySound,'value',1);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',0);
        set(handles.rbScripts,'enable','off');
        
        set(handles.chParam1,'string','one-third-octave-band-spectrogram'); % 440 x 28 x 440
        set(handles.chParam1,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam2,'string','one-third-octave-band-spectrum'); % 1 x 28
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','specific-loudness'); % specific-loudness 1 x 240
        set(handles.chParam3,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
        
    case 12
        
        set(handles.txtAnalysisStart,'Enable','off');
        
        set(handles.rbScripts,'enable','on');
        
        set(handles.rbPsySound,'value',1);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',0);
        
        set(handles.chParam1,'string','loudness');
        set(handles.chParam1,'enable','on');
        
        set(handles.chParam2,'string','main-loudness');
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','specific-loudness');
        set(handles.chParam3,'enable','off');
        set(handles.chParam3,'value',0);
        
        set(handles.chParam4,'string','average-main-loudness');
        set(handles.chParam4,'enable','off');
        set(handles.chParam3,'value',0);
        
        set(handles.chParam5,'string','average-specific-loudness');
        set(handles.chParam5,'enable','on');
        
        set(handles.chParam6,'string','sharpness');
        set(handles.chParam6,'enable','on');
        
        set(handles.chParam7,'string','loudness-fluctuation');
        set(handles.chParam7,'enable','on');
        set(handles.chParam7,'value',0);
    
    case 15
        
        set(handles.txtAnalysisStart,'Enable','on');
        
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'enable','on');
        
        set(handles.chParam1,'string','roughness');
        set(handles.chParam1,'enable','on');
        
        set(handles.chParam2,'string','specific-roughness');
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','average-specific-roughness'); % determined using specific-roughness
        set(handles.chParam3,'enable','on');
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
    
    case 20
        
        set(handles.txtAnalysisStart,'Enable','on');
        
        set(handles.rbPsySound,'enable','off');
        set(handles.rbScripts,'value',1);
        
        set(handles.rbScripts,'enable','on');
        
        set(handles.chParam1,'string','fluctuation-strength');
        set(handles.chParam1,'enable','on');
        
        set(handles.chParam2,'string','specific-fluctuation-strength');
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','average-specific-fluctuation-strength'); % determined using specific-fluctuation-strength
        set(handles.chParam3,'enable','on');
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
        
    case 21
        
        set(handles.txtAnalysisStart,'Enable','off');
        
        set(handles.rbScripts,'value',1);
        set(handles.rbPsySound,'enable','off');
        
        set(handles.chParam1,'string','fundamental-frequency');
        set(handles.chParam1,'enable','on');
        
        set(handles.chParam2,'string','Param2');
        set(handles.chParam2,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','Param3'); 
        set(handles.chParam3,'enable','off');
        
        set(handles.chParam4,'string','Param4');
        set(handles.chParam4,'enable','off');
        set(handles.chParam4,'value',0);
        
        set(handles.chParam5,'string','Param5');
        set(handles.chParam5,'enable','off');
        set(handles.chParam5,'value',0);
        
        set(handles.chParam6,'string','Param6');
        set(handles.chParam6,'enable','off');
        set(handles.chParam6,'value',0);
        
        set(handles.chParam7,'string','Param7');
        set(handles.chParam7,'enable','off');
        set(handles.chParam7,'value',0);
        
    otherwise
        for i = 1:7
            exp1 = sprintf('set(handles.chParam%.0f,''string'',''Param%.0f''); set(handles.chParam%.0f,''enable'',''off''); set(handles.chParam%.0f,''value'',0); ',i,i,i,i);
            eval( exp1 );
        end
        
end

unitgroup_SelectionChangeFcn(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popAnalyser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popAnalyser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtOutputDir_Callback(hObject, eventdata, handles)
% hObject    handle to txtOutputDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOutputDir as text
%        str2double(get(hObject,'String')) returns contents of txtOutputDir as a double


% --- Executes during object creation, after setting all properties.
function txtOutputDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOutputDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtFile1_Callback(hObject, eventdata, handles)
% hObject    handle to txtFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFile1 as text
%        str2double(get(hObject,'String')) returns contents of txtFile1 as a double


% --- Executes during object creation, after setting all properties.
function txtFile1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFile2_Callback(hObject, eventdata, handles)
% hObject    handle to txtFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFile2 as text
%        str2double(get(hObject,'String')) returns contents of txtFile2 as a double


% --- Executes during object creation, after setting all properties.
function txtFile2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in chParam1.
% function chParam1_Callback(hObject, eventdata, handles)
% % hObject    handle to chParam1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of chParam1
% 
% 
% % --- Executes on button press in chParam2.
% function chParam2_Callback(hObject, eventdata, handles)
% % hObject    handle to chParam2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of chParam2
% 
% 
% % --- Executes on button press in chParam3.
% function chParam3_Callback(hObject, eventdata, handles)
% % hObject    handle to chParam3 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of chParam3
% 
% 
% % --- Executes on button press in chParam4.
% function chParam4_Callback(hObject, eventdata, handles)
% % hObject    handle to chParam4 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of chParam4

%% Load data
%   - Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dir_out = get(handles.txtOutputDir,'string');

if length(dir_out) == 0
    
    try
        dir_out = Get_TUe_paths('outputs');
    catch
        warning('Type your output dir in the GUI');
    end
    set(handles.txtOutputDir,'string',dir_out)
    
else
    
    if ~strcmp( dir_out(end), delim )
        dir_out = [dir_out delim];
    end
    
    Mkdir(dir_out);
    
end

filename1 = get(handles.txtFile1,'string');
filename2 = get(handles.txtFile2,'string');

if strcmp(filename1,'')|strcmp(filename2,'')
    try
        filename1 = [dir_out 'tmp-cal' delim 'ref_loud.wav'];
        filename2 = [dir_out 'tmp-cal' delim 'ref_rough.wav'];
    catch
        warning('Enter you wav filenames...')
    end

    set(handles.txtFile1    ,'string',filename1);
    set(handles.txtFile2    ,'string',filename2);
end
set(handles.txtOutputDir,'string',dir_out)

handles.audio.G1 = str2num( get(handles.txtGain1,'string') );
handles.audio.G2 = str2num( get(handles.txtGain2,'string') );

[x1,fs1] = Wavread(filename1);
[x2,fs2] = Wavread(filename2);

if fs1 == fs2
    
    handles.audio.fs = fs1;
    
end

t1 = ( 0:length(x1)-1 )/fs1;
t2 = ( 0:length(x2)-1 )/fs2;

txt2display = sprintf('Length: %.3f [s],\n\t %.0f [samples]\nSample rate: %.0f [Hz]\n',max(t1),length(x1),fs1); 
set( handles.txtFile1info,'string',txt2display);

txt2display2 = sprintf('Length: %.3f [s],\n\t %.0f [samples]\nSample rate: %.0f [Hz]\n',max(t2),length(x2),fs2); 
set( handles.txtFile2info,'string',txt2display2);

ti = str2num( get(handles.txtti,'string') );
tf = str2num( get(handles.txttf,'string') );
toffset = str2num( get(handles.txtXoffset,'string') );

if length(ti)==0 & length(tf)==0
    ti = 1;
    tf = min( length(x1), length(x2) );
    set( handles.txtti,'string',num2str(ti) );
    set( handles.txttf,'string',num2str(tf) );
end

if ti ~= 1 | (tf ~= length(x1) & tf ~= length(x2) )
    handles.audio.bGenerateExcerpt = 1;
    set(handles.txtExcerpt,'visible','on');
else
    handles.audio.bGenerateExcerpt = 0;
    set(handles.txtExcerpt,'visible','off');
end

xliminf = ti/fs1;
xlimsup = tf/fs1;
set(handles.txtti_s,'string',sprintf('%.3f [s]',xliminf));
set(handles.txttf_s,'string',sprintf('%.3f [s], total time of %.3f',xlimsup,(tf-ti)/fs1));

axes(handles.axes1)
plot(t1,From_dB(handles.audio.G1)*x1);
% title( name2figname( filename1 ) )
xlim([xliminf xlimsup])

axes(handles.axes2)
try
    plot(t2(1:end-toffset+1),From_dB(handles.audio.G2)*x2(toffset:end),'r');
catch
    plot(t2,From_dB(handles.audio.G2)*x2,'r');
end
% title( name2figname( filename2 ) )
xlim([xliminf xlimsup])

callevel = str2num( get(handles.txtCalLevel,'string') );

RMS1 = rmsdb(x1(ti        :tf))+handles.audio.G1+callevel+30; % Zwicker's correction
RMS2 = rmsdb(x2(ti+toffset:tf))+handles.audio.G2+callevel+30; % Zwicker's correction

set( handles.txtRMS1,'string',sprintf('RMS, file 1 = %.2f [dB SPL]',RMS1) )
set( handles.txtRMS2,'string',sprintf('RMS, file 2 = %.2f [dB SPL]',RMS2) )

handles.audio.toffset = toffset;

guidata(hObject,handles)

function txtti_Callback(hObject, eventdata, handles)
% hObject    handle to txtti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtti as text
%        str2double(get(hObject,'String')) returns contents of txtti as a double
ti = str2double(get(hObject, 'String'));
if isnan(ti)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new ti value
handles.audio.ti_samples = ti;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function txtti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txttf_Callback(hObject, eventdata, handles)
% hObject    handle to txttf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tf = str2double(get(hObject, 'String'));
if isnan(tf)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new ti value
handles.audio.tf_samples = tf;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function txttf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txttf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtti_s_Callback(hObject, eventdata, handles)
% hObject    handle to txtti_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtti_s as text
%        str2double(get(hObject,'String')) returns contents of txtti_s as a double


% --- Executes during object creation, after setting all properties.
function txtti_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtti_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txttf_s_Callback(hObject, eventdata, handles)
% hObject    handle to txttf_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txttf_s as text
%        str2double(get(hObject,'String')) returns contents of txttf_s as a double


% --- Executes during object creation, after setting all properties.
function txttf_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txttf_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtCalLevel_Callback(hObject, eventdata, handles)
% hObject    handle to txtCalLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCalLevel as text
%        str2double(get(hObject,'String')) returns contents of txtCalLevel as a double


% --- Executes during object creation, after setting all properties.
function txtCalLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCalLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGain1_Callback(hObject, eventdata, handles)
% hObject    handle to txtGain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGain1 as text
%        str2double(get(hObject,'String')) returns contents of txtGain1 as a double


% --- Executes during object creation, after setting all properties.
function txtGain1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGain2_Callback(hObject, eventdata, handles)
% hObject    handle to txtGain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGain2 as text
%        str2double(get(hObject,'String')) returns contents of txtGain2 as a double


% --- Executes during object creation, after setting all properties.
function txtGain2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLabel1_Callback(hObject, eventdata, handles)
% hObject    handle to txtLabel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLabel1 as text
%        str2double(get(hObject,'String')) returns contents of txtLabel1 as a double


% --- Executes during object creation, after setting all properties.
function txtLabel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLabel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLabel2_Callback(hObject, eventdata, handles)
% hObject    handle to txtLabel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLabel2 as text
%        str2double(get(hObject,'String')) returns contents of txtLabel2 as a double


% --- Executes during object creation, after setting all properties.
function txtLabel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLabel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFreqmin_Callback(hObject, eventdata, handles)
% hObject    handle to txtFreqmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFreqmin as text
%        str2double(get(hObject,'String')) returns contents of txtFreqmin as a double


% --- Executes during object creation, after setting all properties.
function txtFreqmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFreqmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFreqmax_Callback(hObject, eventdata, handles)
% hObject    handle to txtFreqmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFreqmax as text
%        str2double(get(hObject,'String')) returns contents of txtFreqmax as a double


% --- Executes during object creation, after setting all properties.
function txtFreqmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFreqmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtXoffset_Callback(hObject, eventdata, handles)
% hObject    handle to txtXoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

toffset = str2double(get(hObject, 'String'));
ti = str2double(get(handles.txtti, 'String'));

if isnan(toffset)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

if ti+toffset < 1
    set(hObject, 'String', 0);
    errordlg('offset - ti has to be greater than 0','Error');
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function txtXoffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtXoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chParam7.
function chParam7_Callback(hObject, eventdata, handles)
% hObject    handle to chParam7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam7



function txtAnalysisStart_Callback(hObject, eventdata, handles)
% hObject    handle to txtAnalysisStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAnalysisStart as text
%        str2double(get(hObject,'String')) returns contents of txtAnalysisStart as a double


% --- Executes during object creation, after setting all properties.
function txtAnalysisStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAnalysisStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtAnalysisEnd_Callback(hObject, eventdata, handles)
% hObject    handle to txtAnalysisEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAnalysisEnd as text
%        str2double(get(hObject,'String')) returns contents of txtAnalysisEnd as a double


% --- Executes during object creation, after setting all properties.
function txtAnalysisEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAnalysisEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
