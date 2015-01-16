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
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PsySoundControl
% Created on: 16/01/2015
% Last modified on: 16/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
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


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
handles.metricdata.volume = volume;
guidata(hObject,handles)

%% Calculations:
% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Needed: 
%       Num2str
%       ef

bSave       = get(handles.bSave,'value');
nAnalyser   = get(handles.popAnalyser,'value');
dir_output  = get(handles.txtOutputDir,'string');

eval( sprintf('options.bDoAnalyser%s=1;',Num2str(nAnalyser,2)) ); % activates selected processor

filename1 = get(handles.txtFile1,'string');
filename2 = get(handles.txtFile2,'string');

options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
options.callevel = 60; % rms 90 dB SPL = 0 dBFS 

options     = Ensure_field(options,'label1','file1');
options     = Ensure_field(options,'label2','file2');

options.bSave = bSave;
options     = Ensure_field(options,'bPlot',1);
options     = Ensure_field(options,'label','');
options     = Ensure_field(options,'SPLrange',[10 70]);
options     = Ensure_field(options,'frange',[50 5000]);

options     = ef(options,'bDoAnalyser08',0);
options     = ef(options,'bDoAnalyser10',0);
options     = ef(options,'bDoAnalyser11',0);
options     = ef(options,'bDoAnalyser12',0);
options     = ef(options,'bDoAnalyser15',0);

options     = Ensure_field(options,'ylim_bExtend',0);
options     = Ensure_field(options,'ylim_bDrawLine',0);

if options.bSave == 1
    options.dest_folder_fig = dir_output;
end
    
h = []; % handles figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' filename1 ' and ' filename2]);
 
ha1 = [];
ha2 = [];    
tmp_h = [];
tmp_h = [];

options = Ensure_field(options,'calfile',[Get_TUe_paths('db_calfiles') 'track_03.wav']);
options = Ensure_field(options,'callevel',70); % 'AMT' reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.bCosineRamp = 0; % Cos ramp not yet applied for Loudness calculations
options.bCosineRampOnset = 0; %ms

options.bPlot = 0;

options.nAnalyser = nAnalyser;

[out_1 tmp_h tmp_ha]   = PsySoundCL(filename1,options);
[out_2 tmp_h tmp_ha]   = PsySoundCL(filename2,options);

param = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser10 == 1
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if options.bDoAnalyser11 == 1
%     param{end+1}        = 'one-third-OB';
%     [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
%     param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
%     
%     if length(stats.idx) == length(stats.t)
%         legend( sprintf('tot = %.2f dB(Z)',dbsum( out_1.DataSpecOneThirdAvg )) ,...
%                 sprintf('tot = %.2f dB(Z)',dbsum( out_2.DataSpecOneThirdAvg )) );
%     else
%         legend( sprintf('tot = %.2f dB(Z)',dbsum( dbmean( out_1.DataSpecOneThird(stats.idx,:) )) ),...
%                 sprintf('tot = %.2f dB(Z)',dbsum( dbmean( out_2.DataSpecOneThird(stats.idx,:) )) ) );
%     end
%     
%     if isfield(options,'SPLrange')
%         set(gca,'YLim',options.SPLrange);
%         plot(options.frange,[65 65],'k'); % horizontal line
%     end
%     
%     if isfield(options,'frange')
%         set(gca,'XLim',options.frange);
%     end
%     hold on
%     
% end

for i = 1:6
    exp1 = sprintf('bPlotParam%.0f = get(handles.chParam%.0f,''value''); labelParam%.0f = get(handles.chParam%.0f,''string'');',i,i,i,i);
    eval( exp1 );
end

  
if bPlotParam1
    % Loudness, Roughness
    param{end+1}        = labelParam1;
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
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
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

end

if bPlotParam3
    param{end+1}        = labelParam3;
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

if bPlotParam4
    param{end+1}        = labelParam4;
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

if bPlotParam5
    param{end+1}        = labelParam5;
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

if bPlotParam6
    param{end+1}        = labelParam6;
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try % Percentiles for Loudness and Sharpness
%     
% output.Lt       = out_1.t;
% output.L_in1    = out_1.DataLoud;
% output.L_in2    = out_2.DataLoud;
% %%%%
% % Percentiles for Loudness
% one_period_s        = diff(options.tanalysis) / (options.N_periods2analyse+1);
% one_period_in_samples = ceil( one_period_s/min(diff(out_1.t)) );
% N_periods = options.N_periods2analyse+1;
% 
% pLmeas  = Get_percentiles_per_period(out_1.DataLoud,one_period_in_samples);
% pLmodel = Get_percentiles_per_period(out_2.DataLoud,one_period_in_samples);
% 
% output.pL_in1   = pLmeas;
% output.pL_in2   = pLmodel;
% 
% %%%%
% % Percentiles for Sharpness
% one_period_in_samples = ceil( one_period_s/min(diff(out_1.t)) );
% yy1 = buffer(out_1.DataSharp,one_period_in_samples,0);
% yy2 = buffer(out_2.DataSharp,one_period_in_samples,0);
% 
% y = yy1(:,1:N_periods);
% p5 = percentile(y,5);
% p50 = percentile(y,50);
% p95 = percentile(y,95);
% 
% p5_mod = percentile(y,5);
% p50_mod = percentile(y,50);
% p95_mod = percentile(y,95);
% 
% catch
%     warning('Percentile calculation not succeeded, maybe not every Analyser is enabled')
% end
% 
% try % Roughness
% 
% output.Rt       = out_1_15.t;
% output.R_in1    = out_1_15.DataRough;
% output.R_in2    = out_2_15.DataRough;
% 
% % Percentiles for Roughness
% 
% one_period_s        = diff(options.tanalysis) / (options.N_periods2analyse+1);
% one_period_in_samples = ceil( one_period_s/min(diff(out_1_15.t)) );
% 
% pRmeas  = Get_percentiles_per_period(out_1_15.DataRough,one_period_in_samples);
% pRmodel = Get_percentiles_per_period(out_2_15.DataRough,one_period_in_samples);
% 
% output.pR_in1   = pRmeas;
% output.pR_in2   = pRmodel;
% 
% catch
%     warning('Percentile calculation not succeeded, maybe not every Analyser is enabled')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output.h = h; % figure handles

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

disp('')

%%
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.rbPsySound)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

%% Initialisation
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

set(handles.unitgroup, 'SelectedObject', handles.rbPsySound);

try
    set(handles.txtOutputDir,'string',Get_TUe_paths('outputs'))
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
        
        set(handles.chParam1,'string','spectrogram'); % 381 x 1024 x 381
        set(handles.chParam1,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam2,'string','average-power-spectrum'); % 1 x 1024
        set(handles.chParam2,'enable','on');
        % set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','spectral-centroid'); % 
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
    
    case 10
        
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
    
    case 11
        
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
        
    case 12
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
    
    case 15
        
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
        
    otherwise
        for i = 1:6
            exp1 = sprintf('set(handles.chParam%.0f,''string'',''Param%.0f''); set(handles.chParam%.0f,''enable'',''off''); set(handles.chParam%.0f,''value'',0); ',i,i,i,i);
            eval( exp1 );
        end
        
end

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


% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    dir_out = Get_TUe_paths('outputs');
catch
    warning('Type your output dir in the GUI');
end

set(handles.txtOutputDir,'string',dir_out)

try
    filename1 = [dir_out 'tmp-cal' delim 'ref_loud.wav'];
    filename2 = [dir_out 'tmp-cal' delim 'ref_rough.wav'];
catch
    warning('Enter you wav filenames...')
end

set(handles.txtOutputDir,'string',dir_out)
set(handles.txtFile1    ,'string',filename1);
set(handles.txtFile2    ,'string',filename2);
