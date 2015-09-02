function varargout = AMTControl(varargin)
% function varargout = AMTControl(varargin)
%
% 1. Description: 
%       AMTCONTROL MATLAB code for AMTControl.fig
%      
%       AMTCONTROL, by itself, creates a new AMTCONTROL or raises the 
%       existing singleton*.
%
%       H = AMTCONTROL returns the handle to a new AMTCONTROL or the handle 
%       to the existing singleton*.
%
%       AMTCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%       function named CALLBACK in AMTCONTROL.M with the given input arguments.
%
%       AMTCONTROL('Property','Value',...) creates a new AMTCONTROL or raises
%       the existing singleton*.  Starting from the left, property value pairs are
%       applied to the GUI before AMTControl_OpeningFcn gets called.  An
%       unrecognised property name or invalid value makes property application
%       stop.  All inputs are passed to AMTControl_OpeningFcn via varargin.
%
%       *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%       instance to run (singleton)".
% 
%       Line    Stage                               Last updated on
%       43      0. Initialisation                   31/07/2015
%       63      1. AMTControl_OpeningFcn            31/07/2015
%       101     3. btnLoad                          07/08/2015
%       338     4. btnGetTemplate                   07/08/2015
%       433     5. btnSimulateAFC                   10/08/2015
%       793     6. btnCalculate                     24/08/2015
%       1101    7. txtXoffset: waveforms shift      31/07/2015
%       1813    8. popAnalyser_Callback             24/08/2015
%       
% TO DO:
%       Dau1997_1Ch
%       str structure to generate automatic MATLAB script
%       Put FontSize as parameter
%       Create menu, with load parameters
% 
% Created on        : 30/07/2015
% Last modified on  : 24/08/2015
% Last used on      : 24/08/2015 % Remember to check compatibility with template_PsySoundCL.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AMTControl_OpeningFcn, ...
                   'gui_OutputFcn',  @AMTControl_OutputFcn, ...
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
% 1. Executes just before AMTControl is made visible.
function AMTControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AMTControl (see VARARGIN)

% Choose default command line output for AMTControl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AMTControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Outputs from this function are returned to the command line.
function varargout = AMTControl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

txtFreqmin_Callback(handles.txtFreqmin,[],handles);
txtFreqmax_Callback(handles.txtFreqmax,[],handles);

dir_out = get(handles.txtOutputDir,'string');

% if output directory is not specified:
if length(dir_out) == 0 
    
    try
        dir_out = Get_TUe_paths('outputs');
    catch
        warning('Type your output dir in the GUI');
    end
    set(handles.txtOutputDir,'string',dir_out)

% if output directory it is specified:
else
    
    % it checks whether last character is separator '\' (win) or '/' (unix):
    if ~strcmp( dir_out(end), delim )
        dir_out = [dir_out delim];
    end
    
    % creates folder in case it does not exist:
    Mkdir(dir_out);
    
end

direx = [Get_TUe_paths('outputs') 'AMTControl-examples' delim];
filename1 = get(handles.txtFile1,'string');
filename2 = get(handles.txtFile2,'string');

if strcmp(filename1,'')|strcmp(filename2,'')
    try
        filename1 = [direx 'dau1996b_expI_noisemasker.wav'];
        filename2 = [direx 'dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav'];
    catch
        warning('Enter you wav filenames...')
    end

    set(handles.txtFile1    ,'string',filename1);
    set(handles.txtFile2    ,'string',filename2);
end
set(handles.txtOutputDir,'string',dir_out)

G1 = str2num( get(handles.txtGain1,'string') );
G2 = str2num( get(handles.txtGain2,'string') );

[insig1,fs]  = Wavread(filename1);
[insig2,fs2] = Wavread(filename2);

% This sould be the normal case:
if fs == fs2
    handles.audio.fs = fs;
end

if strcmp( get(handles.txtti,'String'), '')
    set(handles.txtti,'String','1');
end
if strcmp(get(handles.txttf,'String'),'' )
    set(handles.txttf,'String',num2str(min(length(insig1),length(insig2))));
end

t1 = ( 0:length(insig1)-1 )/fs;
t2 = ( 0:length(insig2)-1 )/fs2;

txt2display = sprintf('Length: %.3f [s],\n\t %.0f [samples]\nSample rate: %.0f [Hz]\n',max(t1),length(insig1),fs); 
set( handles.txtFile1info,'string',txt2display);

txt2display2 = sprintf('Length: %.3f [s],\n\t %.0f [samples]\nSample rate: %.0f [Hz]\n',max(t2),length(insig2),fs2); 
set( handles.txtFile2info,'string',txt2display2);

txtti_Callback(handles.txtti,[],handles);
txttf_Callback(handles.txttf,[],handles);

try
    ti_samples  = handles.audio.ti_samples;
    tf_samples  = handles.audio.tf_samples;
catch
    % Redundant, but first time running needs it (aparently)
    handles.audio.ti_samples = str2num( get(handles.txtti,'String') );
    handles.audio.tf_samples = str2num( get(handles.txttf,'String') );
    ti_samples  = handles.audio.ti_samples;
    tf_samples  = handles.audio.tf_samples;
end
toffset     = str2num( get(handles.txtXoffset,'string') );

if length(ti_samples)==0 & length(tf_samples)==0
    ti_samples = 1;
    tf_samples = min( length(insig1), length(insig2) );
    set( handles.txtti,'string',num2str(ti_samples) );
    set( handles.txttf,'string',num2str(tf_samples) );
    txtti_Callback(handles.txtti,[],handles);
    txttf_Callback(handles.txttf,[],handles);
    handles.audio.ti_samples = ti_samples;
    handles.audio.tf_samples = tf_samples;
end

if ti_samples ~= 1 | (tf_samples ~= length(insig1) & tf_samples ~= length(insig2) )
    handles.audio.bGenerateExcerpt = 1;
    set(handles.txtExcerpt,'visible','on');
else
    handles.audio.bGenerateExcerpt = 0;
    set(handles.txtExcerpt,'visible','off');
end
bGenerateExcerpt = handles.audio.bGenerateExcerpt;

xliminf = ti_samples/fs;
xlimsup = tf_samples/fs;

axes(handles.axes1)
plot(t1,From_dB(G1)*insig1);
% title( name2figname( filename1 ) )
xlim([xliminf xlimsup])
set(gca,'XTickLabel',''); % Time scale will be the same as in File2 (Axes2)

axes(handles.axes2)
try
    plot(t2(1:end-toffset+1),From_dB(G2)*insig2(toffset:end),'r');
catch
    plot(t2,From_dB(G2)*insig2,'r');
end
xlim([xliminf xlimsup])

set(gca,'YTick',[])

lvl_m_30_dBFS = str2num( get(handles.txtCalLevel,'string') );
adjustmentvalue = 70-lvl_m_30_dBFS; % values calibrated to 100 dB RMS = 0 dBFS

insig1_orig = From_dB(adjustmentvalue+G1) * insig1;    
insig2_orig = From_dB(adjustmentvalue+G2) * insig2;

insig1 = insig1_orig(ti_samples        :tf_samples);
insig2 = insig2_orig(ti_samples+toffset:tf_samples);

if bGenerateExcerpt
    Nextra = length(insig1)-(tf_samples - ti_samples)-1;
    if Nextra >= 8192
        Nextra = 8192;
    end
 
    try
        insig1tmp = insig1_orig( ti_samples:tf_samples + Nextra ); % one additional frame
        insig2tmp = insig2_orig( ti_samples + toffset:tf_samples + toffset + Nextra ); % one additional frame
        insig1 = insig1tmp;
        insig2 = insig2tmp;
    catch
        insig1 = insig1_orig( ti_samples:tf_samples ); % one additional frame
        insig2 = insig2_orig( ti_samples + toffset:tf_samples + toffset ); % one additional frame
        warning('using catch...')
    end
    set(handles.txtExcerpt,'visible','on');

    %% Extracts excerpt:
    fname1_excerpt = [Delete_extension(filename1,'wav') '-e.wav']; 
    fname2_excerpt = [Delete_extension(filename2,'wav') '-e.wav']; 

    Wavwrite(insig1,fs,fname1_excerpt);
    Wavwrite(insig2,fs,fname2_excerpt);

    filename1 = fname1_excerpt;
    filename2 = fname2_excerpt;

    handles.audio.filename1 = filename1;
    handles.audio.filename2 = filename2;
   
else
    set(handles.txtExcerpt,'visible','off');

    handles.audio.filename1 = filename1;
    handles.audio.filename2 = filename2;

end

dBFS = lvl_m_30_dBFS+30+adjustmentvalue;
RMS1 = rmsdb(insig1) + dBFS; % Zwicker's correction
RMS2 = rmsdb(insig2) + dBFS; 

thres_silence = 1/3; % one-third of the median value of the envelope 
[xx xx xx RMS1nosil] = Rmssilence(insig1,fs,thres_silence);
[xx xx xx RMS2nosil] = Rmssilence(insig2,fs,thres_silence);

set( handles.txtRMS1,'string',sprintf('RMS1 = %.1f [dB SPL] (%.1f no silent)',RMS1,RMS1nosil + dBFS) );
set( handles.txtRMS2,'string',sprintf('RMS2 = %.1f [dB SPL] (%.1f no silent)',RMS2,RMS2nosil + dBFS) );

handles.audio.toffset = toffset;

bDFT = get(handles.chAddDFT,'value');

if bDFT
    windowtype = 'hanning';
    dBFS = lvl_m_30_dBFS + 30;
    
    K  = length(insig1)/2;
    [xx y1dB f] = freqfft2(insig1,K,fs,windowtype,dBFS);   
    [xx y2dB  ] = freqfft2(insig2,K,fs,windowtype,dBFS);
    
    axes(handles.axesDFT);
    plot(f,y1dB,'b',f,y2dB,'r'); grid on
    min2plot = max( min(min([y1dB y2dB])),0 ); % minimum of 0 dB
    max2plot = max(max([y1dB y2dB]))+5; % maximum exceeded in 5 dB
    ylim([min2plot max2plot])
    xlabel('Frequency [Hz]')
    ylabel('Log-spectrum [dB]')
    
    nTicks = get(gca,'YTick');
    
    if length(nTicks) == 2 % then we add a third Tick
        
        nTicks = [nTicks(1) mean(nTicks) nTicks(2)];
        set(gca,'YTick',nTicks);
        
    end
    
end

% if get(handles.popNavg,'value') ~= 1 % then file 1 should be 'running'
   handles.audio.insig1orig = insig1_orig;     
% end
handles.audio.insig1 = insig1;
handles.audio.insig2 = insig2;
handles.audio.G1 = G1;
handles.audio.G2 = G2;
guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Executes on button press in btnGetTemplate.
function btnGetTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to btnGetTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    insig1 = handles.audio.insig1;
    insig2 = handles.audio.insig2;
    fs = handles.audio.fs;
catch
    error('Load any data before pressing this button...')
end

Gain4supra = il_get_value_numericPop(handles.popGainSupra);

nAnalyser = il_get_nAnalyser(handles.popAnalyser);

if nAnalyser == 100 | nAnalyser == 101 | nAnalyser == 103 | nAnalyser == 104
    
    fc2plot_idx = get(handles.popFc ,'value'); 
    fc2plot_idx2= get(handles.popFc2,'value'); 
    
    fc_ERB      = fc2plot_idx+2;
    Nelements   = length(get(handles.popFc,'String'));
    
    if fc2plot_idx == Nelements
        disp('Choose an fc value to proceed with the analysis...');
        pause();
        fc2plot_idx = get(handles.popFc ,'value');
    end
    
    if fc2plot_idx2~= Nelements
        fc2plot_idx = fc2plot_idx:fc2plot_idx2;
    end
    
    if length(fc2plot_idx) == 1
        fc      = audtofreq(fc_ERB,'erb'); % used as input for single-channel modelling
    else
        fcmin   = audtofreq(min(fc2plot_idx)+2,'erb'); 
        fcmax   = audtofreq(max(fc2plot_idx)+2,'erb');
    end
    mu      = 0;
    
    bSigma = get(handles.chInternalNoise,'value');
    if bSigma
        tmp     = get(handles.popInternalNoise,'string');
        tmp_idx = get(handles.popInternalNoise,'value');
        sigma   = str2num( tmp{tmp_idx} );
    else
        sigma = 0;
    end
    
    if length(fc2plot_idx) == 1
        bSingleChannel = 1;
        bMultiChannel = 0;
        fc2plot_idx_1 = 1;
    else
        bSingleChannel = 0;
        bMultiChannel = 1;
        fc2plot_idx_1 = 1:length(fc2plot_idx);
    end
    
end

template_test = [];
Ntimes = il_get_nAnalyser(handles.popNavg);
if Ntimes == 1
    bDeterministic = 1;
else
    bDeterministic = 0;
end

bUseRamp  = get(handles.chRampMasker,'value');
bUseRampS = get(handles.chRampSignal,'value');

out_1 = [];
out_2 = [];

sigma = 0;
warning('Sigma with temporal assignment')

if bUseRamp;  tmp.masker_ramp_ms = 20; else; tmp.masker_ramp_ms = 0; end % ramp time in ms
if bUseRampS; tmp.signal_ramp_ms = 20; else; tmp.signal_ramp_ms = 0; end % ramp time in ms
        
switch nAnalyser
    case 100
        
        insig2supra = From_dB(Gain4supra) * insig2;
        insig1 = handles.audio.insig1orig;
        
        [out_1Mean out_2Mean fs_intrep] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996',sigma,Ntimes,fc2plot_idx,tmp);
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
    case 101 % Still testing
        
        insig2supra = From_dB(Gain4supra) * insig2;
        insig1 = handles.audio.insig1orig;
        
        [out_1Mean out_2Mean fs_intrep] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1997',sigma,Ntimes,fc2plot_idx,tmp);
                
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
    case 103
        
        insig2supra = From_dB(Gain4supra) * insig2;
        insig1 = handles.audio.insig1orig;
        
        [out_1Mean out_2Mean fs_intrep] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-modfilterbank',sigma,Ntimes,fc2plot_idx,tmp);
                
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
	case 104
        
        insig2supra = From_dB(Gain4supra) * insig2;
        fbstyle = 'lowpass'; 
        
        for i = 1:Ntimes
            
            if bDeterministic == 1
                insig1s0 = handles.audio.insig1;
                insig1s1 = handles.audio.insig1;
            else
                insig1s0 = il_randomise_insig( handles.audio.insig1orig );
                insig1s0 = Do_cos_ramp( insig1s0(1:length(insig2)),fs,tmp.masker_ramp_ms );
                insig1s1 = il_randomise_insig( handles.audio.insig1orig );
                insig1s1 = Do_cos_ramp( insig1s1(1:length(insig2)),fs,tmp.masker_ramp_ms );
            end
            
            if bMultiChannel
                error('Not implemented yet (on 26/08/2015)')
            end
            if bSingleChannel
                [out_1pre , fc, xx, IntRep] = jepsen2008preproc_1Ch(insig1s0  ,fs,fc,fbstyle);
                [out_2pre , fc] = jepsen2008preproc_1Ch(insig1s1 + insig2supra,fs,fc,fbstyle);
                fs_intrep = IntRep.fs_intrep;
            end
            
            if bDeterministic == 1 
                % Deterministic noise
                [out_1 noise] = Add_gaussian_noise_deterministic(out_1pre,mu,sigma); 
                out_2 = out_2pre+noise; % deterministic noise
                
            else
                % 'Running' noise
                tmp = Add_gaussian_noise(out_1pre(:),mu,sigma); % tmp = Add_gaussian_noise(out_1pre(:,fc2plot_idx),mu,sigma);
                out_1 = [out_1 tmp(:)]; 
                tmp = Add_gaussian_noise(out_2pre(:),mu,sigma); % tmp = Add_gaussian_noise(out_2pre(:,fc2plot_idx),mu,sigma);
                out_2 = [out_2 tmp(:)]; % Add internal noise
            end
            
        end
        
        tmp.fs = IntRep.fs_intrep;
        if bDeterministic
            out_1Mean = out_1(:,fc2plot_idx_1);
            out_2Mean = out_2(:,fc2plot_idx_1);
        else
            out_1Mean = mean(out_1,2);
            out_2Mean = mean(out_2,2);
            
            Mtmp = length(fc2plot_idx_1);
            Ntmp = length(out_1Mean)/Mtmp;
            out_1Mean = reshape(out_1Mean,Ntmp,Mtmp);
            out_2Mean = reshape(out_2Mean,Ntmp,Mtmp);
        end
        
        
        %%%
        % if bDeterministic == 1
        %     % deterministic masker
        %     insig1 = handles.audio.insig1;
        % else
        %     insig1 = il_randomise_insig( handles.audio.insig1orig );
        %     % insig1 = insig1(1:length(insig2));
        % end
        % 
        % [out_1Mean out_2Mean fs_intrep] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996',sigma,Ntimes,fc2plot_idx);
        % template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
        %%%
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
end

template = template_test;

t = ( 1:size(out_1Mean,1) )/fs_intrep;

switch nAnalyser
    case 100
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))
        
    case 101
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))
        % figure;
        % opts.bPlot3D = 0;
        % opts.XLabel = 'Time [s]';
        % opts.YLabel = 'Modulation frequency [Hz]';
        % opts.Zlabel = 'Normalised amplitude';
        % t = (1:size(template,1))/fs;
        % warning('variable t only temporally defined')
        % % Mesh(t,mfc,transpose(out_2Mean-out_1Mean),opts)
        % opts.YLim = [min(min(abs(template))) max(max(abs(template)))];
        % Mesh(t,mfc,transpose( template ),opts)
        
    case {103,104}
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))
end

handles.audio.template      = template;
handles.audio.bDeterministic= bDeterministic;
% handles.audio.fc2plot_idx   = fc2plot_idx;
handles.audio.Ntimes        = Ntimes;
handles.audio.fs_intrep     = fs_intrep;

bSave2workspace = get(handles.chSaveTemplate,'value');
if bSave2workspace
    assignin ('base','template',template)
end

guidata(hObject,handles)
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Executes on button press in btnSimulateAFC.
function btnSimulateAFC_Callback(hObject, eventdata, handles)
% hObject    handle to btnSimulateAFC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nAnalyser = il_get_nAnalyser(handles.popAnalyser);

if nAnalyser == 101
    Level_start = 0; % above test-signal level
else
    % Level_start = 10;
    Level_start = il_get_value_numericPop(handles.popGainSupra);
end

Level_step_i    = il_get_value_numericPop(handles.popStepdB);
Reversals_stop  = il_get_value_numericPop(handles.popNreversals);
bHalveStepSize  = get(handles.chStepSizeHalved,'value');

fs = handles.audio.fs;
fs_intrep = handles.audio.fs_intrep;

t1 = str2num( get(handles.txtTimei,'string') ); 
t2 = str2num( get(handles.txtTimef,'string') );
s1 = max(round(t1*fs_intrep), 1); % time 0 is forced to be sample 1
s2 = round(t2*fs_intrep+1);
idxobs = [s1 s2];
Ntimes = handles.audio.Ntimes; %
Nsim = il_get_value_numericPop(handles.popNsim);

try
    insig1 = handles.audio.insig1;
    insig2 = handles.audio.insig2;
    template = handles.audio.template;
    bDeterministic = handles.audio.bDeterministic;
    
    fs = handles.audio.fs;
catch
    error('Load any data and get the template before pressing this button...');
end

bStochastic = ~bDeterministic;

switch nAnalyser
    case {100, 101, 103, 104}
    
        fc2plot_idx = get(handles.popFc ,'value'); 
        fc2plot_idx2= get(handles.popFc2,'value'); 

        fc2plot_idx = fc2plot_idx:fc2plot_idx2;

        if length(fc2plot_idx) == 1
            fc      = audtofreq(fc2plot_idx+2,'erb'); % used as input for single-channel modelling
        else
            fcmin   = audtofreq(min(fc2plot_idx)+2,'erb'); 
            fcmax   = audtofreq(max(fc2plot_idx)+2,'erb');
        end

        mu      = 0;
        
        bSigma = get(handles.chInternalNoise,'value');
        if bSigma
            tmp     = get(handles.popInternalNoise,'string');
            tmp_idx = get(handles.popInternalNoise,'value');
            sigma   = str2num( tmp{tmp_idx} );
        else
            sigma = 0;
        end

        if length(fc2plot_idx) == 1
            bSingleChannel = 1;
            bMultiChannel = 0;
            fc2plot_idx = 1;
        else
            bSingleChannel = 0;
            bMultiChannel = 1;
        end
    
end

Threshold = [];
bUseRamp       = get(handles.chRampMasker,'value');
bUseRampSignal = get(handles.chRampSignal,'value');

if bUseRamp;       rampdn = 20;     end % ramps in ms
if bUseRampSignal; rampsl = rampdn; end

if bUseRampSignal
    fprintf('Introducing %.0f-ms ramps into test signals\n',rampsl);
    insig2 = Do_cos_ramp(insig2,fs,rampsl);
end

for k = 1:Nsim
    
    Level_step = Level_step_i;
    if k > 1
        fprintf('Running simulation: %.0f of %.0f. Last computed threshold=%.2f [dB, relative]\n',k,Nsim,Threshold(k-1));
    else
        fprintf('Running simulation: %.0f of %.0f.\n',k,Nsim);
    end
    
    if bStochastic == 1
        insig1s0 = il_randomise_insig(handles.audio.insig1orig);
        insig1s1 = il_randomise_insig(handles.audio.insig1orig);
        insig1s2 = il_randomise_insig(handles.audio.insig1orig);
    end
    if bDeterministic == 1
        insig1s0 = insig1;
        insig1s1 = insig1;
        insig1s2 = insig1;
    end
    
    insig1s0 = insig1s0( 1:length(insig2) );
    insig1s1 = insig1s1( 1:length(insig2) );
    insig1s2 = insig1s2( 1:length(insig2) );
    if bUseRamp
        if k==1; fprintf('Introducing %.0f-ms ramps into maskers\n',rampdn); end
        insig1s0 = Do_cos_ramp(insig1s0, fs, rampdn);
        insig1s1 = Do_cos_ramp(insig1s1, fs, rampdn);
        insig1s2 = Do_cos_ramp(insig1s2, fs, rampdn);
    end
        
    Level_current   = Level_start;

    Staircase = [];
    Reversals = [];

    dprime     = [];
    dprimetmp = [];
    dprimetmpsign = [];
    changetrend = -1; % decreasing trend is the normal behaviour
    
    nWrong      = 0;
    nCorrect    = 1;
    nReversal   = 0;
    bSucceeded  = 1; % not failed
    
    while (nReversal < Reversals_stop) & (bSucceeded ==  1) % up to line 765

        Gain2apply = From_dB(Level_current);
        insig2test = Gain2apply * insig2;
        interval1 = insig1s0; % Only noise
        interval2 = insig1s1 + insig2test; % Current signal
        interval1s2 = insig1s2;
            
        switch nAnalyser
            case {100, 101, 103, 104}

                if nAnalyser == 100;
                    if bMultiChannel
                        [out_interval1 , fc] = dau1996preproc(interval1,fs);
                        [out_interval2 , fc] = dau1996preproc(interval2,fs);
                        [out_interval1s2, fc] = dau1996preproc(interval1s2,fs);
                        
                        out_interval1 = out_interval1(:,fc2plot_idx);
                        out_interval2 = out_interval2(:,fc2plot_idx);
                        out_interval1s2 = out_interval1s2(:,fc2plot_idx);
                    end
                    if bSingleChannel
                        [out_interval1] = dau1996preproc_1Ch(interval1,fs,fc);
                        [out_interval2] = dau1996preproc_1Ch(interval2,fs,fc);
                        [out_interval1s2] = dau1996preproc_1Ch(interval1s2,fs,fc);
                    end
                    
                elseif nAnalyser == 101
                    [out_interval1 , fc] = dau1997preproc_1Ch(interval1,fs,fc);
                    [out_interval2 , fc] = dau1997preproc_1Ch(interval2,fs,fc);
                    
                elseif nAnalyser == 103
                    if bMultiChannel
                        [out_interval1] = jepsen2008preproc_multi(interval1,fs,fcmin,fcmax,'resample_intrep');
                        [out_interval2] = jepsen2008preproc_multi(interval2,fs,fcmin,fcmax,'resample_intrep');
                        [out_interval1s2] = jepsen2008preproc_multi(interval1s2,fs,fcmin,fcmax,'resample_intrep');
                        
                        out_interval1   = il_pool_in_one_column(out_interval1);
                        out_interval2   = il_pool_in_one_column(out_interval2);
                        out_interval1s2 = il_pool_in_one_column(out_interval1s2);
                    end
                    if bSingleChannel
                        [out_interval1] = jepsen2008preproc_1Ch(interval1,fs,fc);
                        [out_interval2] = jepsen2008preproc_1Ch(interval2,fs,fc);
                        [out_interval1s2] = jepsen2008preproc_1Ch(interval1s2,fs,fc);
                        
                        out_interval1   = il_pool_in_one_column(out_interval1);
                        out_interval2   = il_pool_in_one_column(out_interval2);
                        out_interval1s2 = il_pool_in_one_column(out_interval1s2);
                    end
                    
                elseif nAnalyser == 104
                    if bMultiChannel
                        error('Not implemented yet (on 26/08/2015)');
                    end
                    if bSingleChannel
                        [out_interval1] = jepsen2008preproc_1Ch(interval1,fs,fc,'lowpass');
                        [out_interval2] = jepsen2008preproc_1Ch(interval2,fs,fc,'lowpass');
                        [out_interval1s2] = jepsen2008preproc_1Ch(interval1s2,fs,fc,'lowpass');
                    end
                end
                
                out_interval1 = Add_gaussian_noise(out_interval1,mu,sigma); % Add internal noise
                out_interval2 = Add_gaussian_noise(out_interval2,mu,sigma); % Add internal noise
                out_interval1s2 = Add_gaussian_noise(out_interval1s2,mu,sigma);
                
                [a b] = size(template);
                % one audio-frequency band but all the modulation filterbanks:
                try
                    sigint1 = reshape(out_interval1,a,b);
                    sigint2 = reshape(out_interval2,a,b);
                    sigint1s2 = reshape(out_interval1s2,a,b);
                catch
                    error('Continue here')
                    sigint1 = out_interval1(idxobs(1):idxobs(2),:);
                    sigint2 = out_interval2(idxobs(1):idxobs(2),:);
                    sigint1s2 = out_interval1s2(idxobs(1):idxobs(2),:);
                end
        end
        
        stats1.mean = mean(sigint1(:));
        stats2.mean = mean(sigint2(:));
        dprime(end+1,1) = (stats2.mean - stats1.mean )/sigma;

        try
            tmp_mue1 = optimaldetector(  sigint1,template);
            tmp_mue2 = optimaldetector(sigint1s2,template);
        catch
            error('continue here')
            tmp_mue1 = optimaldetector(  sigint1,template(idxobs(1):idxobs(2),:));
            tmp_mue2 = optimaldetector(sigint1s2,template(idxobs(1):idxobs(2),:));
        end
        decision(1) = max([tmp_mue1 tmp_mue2]);
        
        decsign = sign([tmp_mue1 tmp_mue2]);
        [dec2m idxtmp] = max([abs(tmp_mue1) abs(tmp_mue2)]); 
        % dec2m = dec2m*decsign(idxtmp);
        
        try
            decision(2) = optimaldetector(sigint2,template);
        catch
            error('continue here')
            decision(2) = optimaldetector(sigint2,template(idxobs(1):idxobs(2),:));
        end
        
        finaldecision = ( decision(2)-decision(1) )/sigma;
        
        fdec = (decision(2)-dec2m)/sigma;
        dprimetmp(end+1) = finaldecision;
                
        if abs( dprime(end) ) < 1.26
            disp('');
        end

        Staircase = [Staircase; Level_current];
        
        % if dprime(end) >= 1.26
        if (abs(dprimetmp(end)) >= 1.26) % & (bHasChanged == 0)
            if nCorrect == 0

                nReversal = nReversal + 1;
                Reversals = [Reversals; Level_current];

                if mod(nReversal,2) == 0 & bHalveStepSize
                    Level_step = max( Level_step/2, 1 );
                end

            end
            nCorrect = 1;
            
            Level_current = Level_current - Level_step; % we make it more difficult
            
        else % masker is chosen
            nWrong = nWrong+1;
            if nWrong == 2

                nWrong = 0;
                nReversal = nReversal + 1;
                Reversals = [Reversals; Level_current];
                if mod(nReversal,2) == 0 & bHalveStepSize
                    Level_step = max( Level_step/2, 1 );
                end
                Level_current = Level_current + Level_step; % we make it easier after two mistakes
                nCorrect = 0;
            end
        end

        if (Level_current < -120)
            bSucceeded = 0;
        else
            bSucceeded = 1;
        end
        
        if mod(length(Staircase),60) == 0
            fprintf('Number of simlulated trials for this run so far = %.0f\n',length(Staircase));
        end

    end
    
    if bSucceeded
        Threshold(k) = median(Reversals(end-6+1:end,:));
    else
        Threshold(k) = NaN;
    end
    
    if bDeterministic
        figure; 
        plot(Staircase,'o');
        xlabel('Presentation order')
        ylabel('Relative level')
        title(sprintf('Reversals median = %.2f [dB]',Threshold(k)));
        grid on
    end
    
    if bStochastic & k == Nsim
        figure;
        plot(Threshold,'o');
        xlabel('Running noise nr.')
        ylabel('Level at threshold (relative level)')
        title(sprintf('Average threshold median (L25,L75) = %.2f dB (%.2f, %.2f)',prctile(Threshold,50),prctile(Threshold,25),prctile(Threshold,75)));
        grid on
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Calculations: 
%   Executes on button press in calculate.
function btnCalculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Needed: 
%       Num2str
%       ef

bUsePsySound = get(handles.rbPsySound,'value');
options.bUsePsySound = bUsePsySound;

bSave = get(handles.bSave,'value');

nAnalyser = il_get_nAnalyser(handles.popAnalyser);

dir_output  = get(handles.txtOutputDir,'string');

nSkipStart  = str2num( get(handles.txtAnalysisStart,'string') );
nSkipEnd    = str2num( get(handles.txtAnalysisEnd  ,'string') );

if nAnalyser == 15 | bUsePsySound
    HopSize = str2num( get(handles.txtOverlap,'string') );
    CParams.HopSize = HopSize;
end

if nAnalyser == 100 | nAnalyser == 101
    
    fc2plot_idx = get(handles.popFc ,'value'); 
    fc_ERB = fc2plot_idx+2;
    fc     = audtofreq(fc_ERB,'erb');
    mu      = 0;
    tmp     = get(handles.popInternalNoise,'string');
    tmp_idx = get(handles.popInternalNoise,'value');
    sigma   = str2num( tmp{tmp_idx} );
    
end

options.bLogScale = get(handles.cbLogAxis,'value'); % to be used in PsySound_Figures

fs         = handles.audio.fs;
sample_inf = handles.audio.ti_samples;
sample_sup = handles.audio.tf_samples;

toffset = handles.audio.toffset; % time offset of audio 2 in relation to audio 1

options.bGenerateExcerpt= handles.audio.bGenerateExcerpt;
bGenerateExcerpt        = handles.audio.bGenerateExcerpt;

if bGenerateExcerpt
    options.tanalysis = [sample_inf-1 sample_sup-1]/fs;
    set( handles.txtExcerpt,'visible','on');
 else
    set( handles.txtExcerpt,'visible','off');
end

eval( sprintf('options.bDoAnalyser%s=1;',Num2str(nAnalyser,2)) ); % activates selected processor

filename1 = handles.audio.filename1; % get(handles.txtFile1,'string');
filename2 = handles.audio.filename2; % get(handles.txtFile2,'string');

options.nAnalyser   = nAnalyser;
options.bSave       = bSave;

if bSave == 1
    options.dest_folder_fig = dir_output;
end
    
% Plot options:
options.label1 = get(handles.txtLabel1,'string');
options.label2 = get(handles.txtLabel2,'string');

% options     = Ensure_field(options,'bPlot',1);
options          = Ensure_field(options,'label','');
options.SPLrange = [str2num(get(handles.txtlvlmin ,'string'))  str2num(get(handles.txtlvlmax ,'string'))];
options.frange   = [str2num(get(handles.txtFreqmin,'string'))  str2num(get(handles.txtFreqmax,'string'))];
options.zrange   = [str2num(get(handles.txtZmin   ,'string'))  str2num(get(handles.txtZmax   ,'string'))];
options.trange   = [str2num(get(handles.txtTimei  ,'string'))  str2num(get(handles.txtTimef  ,'string'))];
options          = Ensure_field(options,'ylim_bExtend',0);
options          = Ensure_field(options,'ylim_bDrawLine',0);
options.bLoudnessContrained = get(handles.chAvgLoudnessLimits,'Value'); % only validated for Analyser 12
options.zlim4assessment = options.zrange; 

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
    options.nSkipEnd   = nSkipEnd;
    [out_1 tmp_h tmp_ha]   = PsySoundCL(filename1,options,CParams);
    
    options.callevel = callevel + handles.audio.G2;
    options.nSkipStart = nSkipStart;
    options.nSkipEnd   = nSkipEnd;
    [out_2 tmp_h tmp_ha]   = PsySoundCL(filename2,options,CParams);
    
else % if NOT PsySound
    
    callevel = str2num( get(handles.txtCalLevel,'string') ); % rms 90 dB SPL = 0 dBFS 
    warning('callevel not used at all for m-file scripts...');
    
    insig1 = handles.audio.insig1; % [insig1 fs1] = Wavread(filename1);
    insig2 = handles.audio.insig2; % [insig2 fs2] = Wavread(filename2);
    
    lvl_m_30_dBFS = str2num( get(handles.txtCalLevel,'string') );
    calvalue = lvl_m_30_dBFS-70; 
    
    switch nAnalyser
        
        case 1
            
            windowtype = 'hanning';
            dBFS = lvl_m_30_dBFS + 30;

            K  = length(insig1)/2;
            [xx y1dB f] = freqfft2(insig1,K,fs,windowtype,dBFS);   
            [xx y2dB  ] = freqfft2(insig2,K,fs,windowtype,dBFS);

            out_1.t = ( 1:length(insig1) )/fs;
            out_1.f = f;
            out_1.Data1 = y1dB;
            
            nParam = 2;
            out_1.name{nParam} = 'Log-spectrum';
            out_1.param{nParam} = strrep( lower( out_1.name{nParam} ),' ','-');
            
            out_2.t = ( 1:length(insig2) )/fs;
            out_2.f = f;
            out_2.Data1 = y2dB;
                        
            nParam = 2;
            out_2.name{nParam} = 'Log-spectrum';
            out_2.param{nParam} = strrep( lower( out_2.name{nParam} ),' ','-');


        case 12 % Loudness
            
            % Only loudness fluctuation:
            dBFS = lvl_m_30_dBFS + 30 - calvalue;
            
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
            opts.nSkipEnd   = nSkipEnd;
            [xx xx out_1] = Roughness_offline(insig1,fs,N,opts,CParams,0);
            [xx xx out_2] = Roughness_offline(insig2,fs,N,opts,CParams,0);
            
        case 20 % Fluctuation strength, see also r20141126_fluctuation
            
            N = 44100*4; % 8192*4;
            opts.nSkipStart = nSkipStart;
            opts.nSkipEnd   = nSkipEnd;
            warning('Fluctuation strength: temporal value...')
            [xx out_1] = FluctuationStrength_offline_debug(insig1(1:N),fs,N,0);
            [xx out_2] = FluctuationStrength_offline_debug(insig2(1:N),fs,N,0);
            
        case 100
            
            bPlotParams = il_get_plotParams(handles);
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_1 , fc] = dau1996preproc(insig1,fs);
                [out_1 xx] = Add_gaussian_noise(out_1,mu,sigma); % Add internal noise
                
            else
                
                % [out_1 , fc, extraouts] = dau1996preproc_1Ch(insig1,fs,fc);
                % if bPlotParams(2) == 1
                %     filteroutsig = extraouts.out_filterbank;
                % 
                %     f1 = sprintf('%sAMTControl-Examples%sfile1-audfilter-fc-%.0f-Hz.wav',Get_TUe_paths('outputs'),delim,fc);
                %     Wavwrite(filteroutsig,fs,f1);
                % end
                [out_1 , fcs, extraouts] = dau1996preproc(insig1,fs); % [out_2 , fc, extraouts] = dau1996preproc_1Ch(insig2,fs,fc);
                idx_tmp = find(fc >= fcs); idx_tmp = idx_tmp(end);
                out_1 = out_1(:,idx_tmp);
                
                filteroutsig = extraouts.out_filterbank;
                dirtmp = sprintf('%sAMTControl-Examples%sGammatone-out%s',Get_TUe_paths('outputs'),delim,delim);
                Mkdir(dirtmp);
                for i = 1:length(fcs)
                    f1 = sprintf('%sfile1-audfilter-fc-%.0f-Hz.wav',dirtmp,fcs(i));
                    Wavwrite(filteroutsig(:,i),fs,f1);
                end
            end
            
            out_1.out   = out_1;
            out_1.fs    = fs;
            out_1.fc    = fc;
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_2 , fc] = dau1996preproc(insig2,fs);
                out_2 = Add_gaussian_noise(out_2,mu,sigma); % Add internal noise
                
            else
                
                [out_2 , fcs, extraouts] = dau1996preproc(insig2,fs);
                idx_tmp = find(fc >= fcs); idx_tmp = idx_tmp(end);
                out_2 = out_2(:,idx_tmp);
                
                filteroutsig = extraouts.out_filterbank;
                dirtmp = sprintf('%sAMTControl-Examples%sGammatone-out%s',Get_TUe_paths('outputs'),delim,delim);
                for i = 1:length(fcs)
                    f1 = sprintf('%sfile2-audfilter-fc-%.0f-Hz.wav',dirtmp,fcs(i));
                    Wavwrite(filteroutsig(:,i),fs,f1);
                end
                
                % if bPlotParams(2) == 1
                %     filteroutsig = extraouts.out_filterbank;
                % 
                %     f1 = sprintf('%sAMTControl-Examples%sfile2-audfilter-fc-%.0f-Hz.wav',Get_TUe_paths('outputs'),delim,fc);
                %     Wavwrite(filteroutsig,fs,f1);
                % end
                
            end
            
            out_2.out   = out_2;
            out_2.fs    = fs;
            out_2.fc    = fc;
            
        case 101
            
            bPlotParams = il_get_plotParams(handles);
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_1 , fc ,fcm] = dau1997preproc(insig1,fs);
                out_1 = Add_gaussian_noise(out_1{fc2plot_idx},mu,sigma); % Add internal noise
                
            else
                
                [out_1 , fc ,fcm, extraouts] = dau1997preproc_1Ch(insig1,fs,fc);
                if bPlotParams(2) == 1
                    filteroutsig = extraouts.out01_filterbank;
                    
                    f1 = sprintf('%sAMTControl-Examples%sfile1-audfilter-fc-%.0f-Hz.wav',Get_TUe_paths('outputs'),delim,fc);
                    Wavwrite(filteroutsig,fs,f1);
                end
            end
            
            out_1.out   = out_1;
            
            out_1.fs    = fs;
            out_1.fc    = fc;
            out_1.fcm   = fcm;
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_2 , fc, fcm] = dau1997preproc(insig2,fs);
                out_2 = Add_gaussian_noise(out_2{fc2plot_idx},mu,sigma); % Add internal noise
                
            else
                
                [out_2 , fc, fcm, extraouts] = dau1997preproc_1Ch(insig2,fs,fc);
                if bPlotParams(2) == 1
                    filteroutsig = extraouts.out01_filterbank;
                    
                    f1 = sprintf('%sAMTControl-Examples%sfile2-audfilter-fc-%.0f-Hz.wav',Get_TUe_paths('outputs'),delim,fc);
                    Wavwrite(filteroutsig,fs,f1);
                end
                
            end
            out_2.out   = out_2;
            out_2.fs    = fs;
            out_2.fc    = fc;
            out_2.fcm   = fcm;
            
        case 104
            
            bPlotParams = il_get_plotParams(handles);
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_1, fc, fcm] = jepsen2008preproc(insig1,fs,'lowpass');
                [out_1 xx] = Add_gaussian_noise(out_1,mu,sigma); % Add internal noise
                
            else
                
                [out_1, fc, fcm, extraouts] = jepsen2008preproc(insig1,fs,'lowpass');
                if bPlotParams(2) == 1
                    filteroutsig = extraouts.out_filterbank;
                    
                    dirtmp = sprintf('%sAMTControl-Examples%sDRNL-out%s',Get_TUe_paths('outputs'),delim,delim);
                    Mkdir(dirtmp);
                    for i = 1:length(fc)
                        f1 = sprintf('%sfile1-audfilter-fc-%.0f-Hz.wav',dirtmp,fc(i));
                        Wavwrite(filteroutsig(:,i),fs,f1);
                    end
                end
            end
            
            out_1.out   = out_1;
            out_1.fs    = fs;
            out_1.fc    = fc;
            out_1.fcm   = fcm;
            
            if sum( bPlotParams(2:4) ) == 0 % we need only the final outout of the peripheral processing
                
                [out_2, fc, fcm] = jepsen2008preproc(insig2,fs,'lowpass');
                out_2 = Add_gaussian_noise(out_2,mu,sigma); % Add internal noise
                
            else
                
                [out_2, fc, fcm, extraouts] = jepsen2008preproc(insig2,fs,'lowpass');
                if bPlotParams(2) == 1
                    filteroutsig = extraouts.out_filterbank;
                    
                    dirtmp = sprintf('%sAMTControl-Examples%sDRNL-out%s',Get_TUe_paths('outputs'),delim,delim);
                    for i = 1:length(fc)
                        f1 = sprintf('%sfile2-audfilter-fc-%.0f-Hz.wav',dirtmp,fc(i));
                        Wavwrite(filteroutsig(:,i),fs,f1);
                    end
                end
            end
            
            out_2.out   = out_2;
            out_2.fs    = fs;
            out_2.fc    = fc;
            out_2.fcm   = fcm;
            
    end
    
end

param   = [];
h       = []; % handles figures
ha      = [];

for i = 1:7
    exp1 = sprintf('bPlotParam%.0f = get(handles.chParam%.0f,''value''); labelParam%.0f = get(handles.chParam%.0f,''string'');',i,i,i,i);
    eval( exp1 );
    
    % To generate automatically the script:
    exp1 = sprintf('str.bPlotParam%.0f = get(handles.chParam%.0f,''value''); str.labelParam%.0f = get(handles.chParam%.0f,''string'');',i,i,i,i);
    eval( exp1 );
end

%% Plots
    
switch nAnalyser
    case 100
        
        if bPlotParams(1) == 1
            optsTmp = [];
            optsTmp.bSave  = 0;
            optsTmp.Title  = sprintf('%s',options.label1); 
            h(end+1) = il_Plot_dau1996(out_1,optsTmp);
            param{end+1}   = sprintf('%s-analyser-%s-file1','MU',Num2str(options.nAnalyser));
                    
            optsTmp = [];
            optsTmp.bSave  = 0;
            optsTmp.Title  = sprintf('%s',options.label2); 
            h(end+1) = il_Plot_dau1996(out_2,optsTmp);
            param{end+1}   = sprintf('%s-analyser-%s-file2','MU',Num2str(options.nAnalyser));
        end
        
    case 101
        
        if length(fc) > 1
            fc2plot = fc(fc2plot_idx);
        else
            fc2plot = fc;
        end
        
        if bPlotParams(1) == 1
            optsTmp = [];
            optsTmp.bSave  = 0; % is going to be saved by this GUI, later
            optsTmp.Title  = sprintf('%s, fc=%.0f [Hz]',options.label1,fc2plot); 
            % optsTmp.YLim   = [-105 400];
            h(end+1) = il_Plot_dau1997(out_1,optsTmp);
            param{end+1}   = sprintf('%s-analyser-%s-file1','MU',Num2str(options.nAnalyser));

            optsTmp = [];
            optsTmp.bSave  = 0;
            optsTmp.Title  = sprintf('%s, fc=%.0f [Hz]',options.label2,fc2plot); 
            % optsTmp.YLim   = [-105 400];
            h(end+1) = il_Plot_dau1997(out_2,optsTmp);
            param{end+1}   = sprintf('%s-analyser-%s-file2','MU',Num2str(options.nAnalyser));
        end
        
    % case 1
    % 
    %     if bUsePsySound == 0
    % 
    %         if bPlotParam2 == 1
    %             figure;
    %             plot(f,y1dB,'b',f,y2dB,'r'); grid on
    %             min2plot = max( min(min([y1dB y2dB])),0 ); % minimum of 0 dB
    %             max2plot = max(max([y1dB y2dB]))+5; % maximum exceeded in 5 dB
    %             ylim([min2plot max2plot])
    %             xlabel('Frequency [Hz]')
    %             ylabel('Log-spectrum [dB]')
    % 
    %             legend(options.label1,options.label2);
    % 
    %         end
    % 
    %     end
        
    otherwise
        
        bPercentiles = get(handles.chPercentiles,'value'); % only for loudness

        if bPercentiles & nAnalyser == 12
            param{end+1} = 'loudness-percentiles';
            [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1,out_2,options);
            param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
        end

        if bPlotParam1
            % Loudness, Roughness
            param{end+1} = labelParam1; % to be used in PsySoundCL_Figures
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
            % Log-spectrum, Specific loudness, roughness
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
        param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

end

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
            
            options.tanalysis = [sample_inf-1 sample_sup-1]/fs;
            
            fndiary = [options.dest_folder_fig 'fig-log-analyser-' num2str(nAnalyser) '.txt'];
            diary(fndiary)

            p = Get_date;
            str.functionname = sprintf('PS_CL_%s_%s',Num2str(nAnalyser),p.date4files); 
            str.ifile = [Get_TUe_paths('MATLAB') 'template_PsySoundCL.m'];
            str.bSave = bSave;
            str.f1 = filename1;
            str.f2 = filename2;
            str.ofile = [Get_TUe_paths('MATLAB') str.functionname '.m'];
            str.nSkipStart = nSkipStart;
            str.nSkipEnd = nSkipEnd;
            str.bUsePsySound = bUsePsySound;
            str.fs = fs;
            str.sample_inf = sample_inf;
            str.sample_sup = sample_sup;
            str.HopSize = HopSize;
            str.toffset = toffset;
            str.bGenerateExcerpt = options.bGenerateExcerpt;
            str.tanalysis = options.tanalysis;
            str.nAnalyser = nAnalyser;
            str.label1 = options.label1;
            str.label2 = options.label2;
            str.SPLrange = options.SPLrange;
            str.frange = options.frange;
            str.zrange = options.zrange;
            str.trange = options.trange;
            str.bLoudnessContrained = options.bLoudnessContrained;
            str.G1 = handles.audio.G1;
            str.G2 = handles.audio.G2;
            str.callevel = callevel;
            try
                str.bPercentiles = bPercentiles;
            end

            o2write = readfile_replace(str.ifile,str);

            ofile = str.ofile;

            fid=fopen(ofile, 'w'); 
            fwrite(fid, o2write); 
            fclose(fid);
        
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf('Date/time of processing: %s\n', Get_date_ddmmyyyy(1));
            fprintf('Output directory: %s\n'    ,options.dest_folder_fig);
            fprintf('Level ref. tone: %.1f dB\n',callevel);
            fprintf('File name 1: %s (gain = %.2f dB)\n',filename1,handles.audio.G1);
            fprintf('File name 2: %s (gain = %.2f dB)\n',filename2,handles.audio.G2);
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
        
    end
    
else
    
    disp('Figures are NOT going to be stored... Set options.bSave to 1 and re-run the scripts in case you want to save the figures')
    pause(1)
    
end

if ~bUsePsySound
    
    try
        delete( fname1_excerpt );
        delete( fname2_excerpt );
        disp('...deleting temporal audio file...');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in cbLogAxis.
function cbLogAxis_Callback(hObject, eventdata, handles)
% hObject    handle to cbLogAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbLogAxis



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


% --- Executes on selection change in txtWindowsSize.
function txtWindowsSize_Callback(hObject, eventdata, handles)
% hObject    handle to txtWindowsSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns txtWindowsSize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from txtWindowsSize


% --- Executes during object creation, after setting all properties.
function txtWindowsSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWindowsSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZmin_Callback(hObject, eventdata, handles)
% hObject    handle to txtZmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZmin as text
%        str2double(get(hObject,'String')) returns contents of txtZmin as a double


% --- Executes during object creation, after setting all properties.
function txtZmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZmax_Callback(hObject, eventdata, handles)
% hObject    handle to txtZmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZmax as text
%        str2double(get(hObject,'String')) returns contents of txtZmax as a double


% --- Executes during object creation, after setting all properties.
function txtZmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtlvlmin_Callback(hObject, eventdata, handles)
% hObject    handle to txtlvlmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtlvlmin as text
%        str2double(get(hObject,'String')) returns contents of txtlvlmin as a double


% --- Executes during object creation, after setting all properties.
function txtlvlmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtlvlmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtlvlmax_Callback(hObject, eventdata, handles)
% hObject    handle to txtlvlmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtlvlmax as text
%        str2double(get(hObject,'String')) returns contents of txtlvlmax as a double


% --- Executes during object creation, after setting all properties.
function txtlvlmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtlvlmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTimei_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimei as text
%        str2double(get(hObject,'String')) returns contents of txtTimei as a double


% --- Executes during object creation, after setting all properties.
function txtTimei_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTimef_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimef as text
%        str2double(get(hObject,'String')) returns contents of txtTimef as a double


% --- Executes during object creation, after setting all properties.
function txtTimef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txtti_Callback(hObject, eventdata, handles)
% hObject    handle to txtti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ti = str2num( get(handles.txtti,'string') );
tf = str2num( get(handles.txttf,'string') );
if isnan(ti)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

try
    fs = handles.audio.fs;
catch
    fs = 44100;
    disp('Assuming theoretical fs')
end

xliminf = (ti-1)/fs;
xlimsup = (tf-1)/fs;

set(handles.txtTimei,'string',num2str(xliminf));
set(handles.txtTimef,'string',num2str(xlimsup));

set(handles.txtti_s,'string',sprintf('%.3f',xliminf));
set(handles.txttf_s,'string',sprintf('%.3f',xlimsup));
set(handles.txtAnalysisTime,'string',sprintf('%.3f',(tf-ti)/fs));

% Save the new ti value
handles.audio.ti_samples = ti;
handles.audio.tf_samples = tf;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txttf_Callback(hObject, eventdata, handles)
% hObject    handle to txttf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ti = str2num( get(handles.txtti,'string') );
tf = str2num( get(handles.txttf,'string') );

if isnan(tf)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

try
    fs = handles.audio.fs;
catch
    fs = 44100;
    disp('Assuming theoretical fs')
end

xliminf = (ti-1)/fs;
xlimsup = (tf-1)/fs;

set(handles.txtTimei,'string',num2str(xliminf));
set(handles.txtTimef,'string',num2str(xlimsup));

set(handles.txtti_s,'string',sprintf('%.3f',xliminf));
set(handles.txttf_s,'string',sprintf('%.3f',xlimsup));
set(handles.txtAnalysisTime,'string',sprintf('%.3f',(tf-ti)/fs));

% Save the new ti value
handles.audio.ti_samples = ti;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X. txtXoffset_Callback: waveforms shift
function txtXoffset_Callback(hObject, eventdata, handles)
%   - hObject    handle to txtXoffset (see GCBO)
%   - eventdata  reserved - to be defined in a future version of MATLAB
%   - handles    structure with handles and user data (see GUIDATA)

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

function txtAnalysisTime_Callback(hObject, eventdata, handles)
% hObject    handle to txtAnalysisTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAnalysisTime as text
%        str2double(get(hObject,'String')) returns contents of txtAnalysisTime as a double


% --- Executes during object creation, after setting all properties.
function txtAnalysisTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAnalysisTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chPercentiles.
function chPercentiles_Callback(hObject, eventdata, handles)
% hObject    handle to chPercentiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chPercentiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Executes on selection change in popAnalyser.
function popAnalyser_Callback(hObject, eventdata, handles)
%   - hObject    handle to popAnalyser (see GCBO)
%   - eventdata  reserved - to be defined in a future version of MATLAB
%   - handles    structure with handles and user data (see GUIDATA)

% nAnalyser = get(handles.popAnalyser,'value');
nAnalyser = il_get_nAnalyser(handles.popAnalyser);

switch nAnalyser

    case 1
        
        set(handles.txtAnalysisStart,'Enable','off');
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.popInternalNoise,'enable','off');
        set(handles.popFc,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','off');
        
        set(handles.rbPsySound,'value',0);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',1);
        set(handles.rbScripts,'enable','on');
        
        set(handles.chParam1,'string','spectrogram'); % 381 x 1024 x 381
        set(handles.chParam1,'enable','on');
        set(handles.chParam1,'value',0);
        
        set(handles.chParam2,'string','log-spectrum'); % 1 x 1024
        set(handles.chParam2,'enable','on');
        set(handles.chParam2,'value',1);
        
        set(handles.chParam3,'string','spectral-centroid'); % 
        set(handles.chParam3,'enable','off');
        set(handles.chParam3,'value',0);
        
        set(handles.chParam4,'string','average-power-spectrum');
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
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.popInternalNoise,'enable','off');
        set(handles.popFc,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','off');
        
        set(handles.rbPsySound,'value',1);
        set(handles.rbPsySound,'enable','on');
        set(handles.rbScripts,'value',0);
        set(handles.rbScripts,'enable','off');
        
        set(handles.chParam1,'string','one-third-octave-band-spectrogram'); % 440 x 28 x 440
        set(handles.chParam1,'enable','off');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam2,'string','one-third-octave-band-spectrum'); % 1 x 28
        set(handles.chParam2,'enable','on');
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
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.popInternalNoise,'enable','off');
        set(handles.popFc,'enable','off');
        set(handles.chPercentiles,'enable','on');
        set(handles.chInternalNoise,'enable','off');
        
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
    
    case 15 % Roughness
        
        set(handles.txtAnalysisStart,'Enable','on');
        set(handles.txtAnalysisEnd  ,'Enable','on');
        set(handles.popInternalNoise,'enable','off');
        set(handles.popFc,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','off');
        
        set(handles.rbPsySound,'enable','off');
        set(handles.rbPsySound,'value',0);
        set(handles.rbScripts,'enable','on');
        set(handles.rbScripts,'value',1);
        
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
    
        % Customable params:
        set(handles.txtOverlap,'string',num2str(4096));
        
    case 20 % Fluctuation strength
        
        set(handles.txtAnalysisStart,'Enable','on');
        set(handles.txtAnalysisEnd  ,'Enable','on');
        set(handles.popInternalNoise,'enable','off');
        set(handles.popFc,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','off');
        
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
        
    case 100 % dau1996
        
        set(handles.txtCalLevel,'string',70);
        set(handles.txtAnalysisStart,'Enable','off');
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.txtOverlap,'enable','off');
        
        set(handles.chAvgLoudnessLimits,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','on');
        set(handles.popInternalNoise,'enable','on');
        set(handles.popFc,'enable','on');
        set(handles.rbScripts,'value',1);
        set(handles.rbPsySound,'enable','off');
                
        set(handles.chParam1,'string','Low-passed filters');
        set(handles.chParam1,'enable','on');
        set(handles.chParam1,'value',1);
        
        set(handles.chParam2,'string','Auditory filterbank');
        set(handles.chParam2,'enable','on');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','Hair-cell envelope'); 
        set(handles.chParam3,'enable','off');
        
        set(handles.chParam4,'string','Adaptation loops');
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

    case 101 % dau1997
        
        set(handles.txtCalLevel,'string',70);
        set(handles.txtAnalysisStart,'Enable','off');
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.txtOverlap,'enable','off');
        
        set(handles.popInternalNoise,'enable','on');
        set(handles.popFc,'enable','on');
        set(handles.chAvgLoudnessLimits,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','on');
        set(handles.rbScripts,'value',1);
        set(handles.rbPsySound,'enable','off');
        
        set(handles.chParam1,'string','Modulation filters');
        set(handles.chParam1,'enable','on');
        set(handles.chParam1,'value',1);
        
        set(handles.chParam2,'string','Auditory filterbank');
        set(handles.chParam2,'enable','on');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','Hair-cell envelope'); 
        set(handles.chParam3,'enable','off');
        
        set(handles.chParam4,'string','Adaptation loops');
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
    
    case 104 % jepsen2008 (LP modulation filter)
        
        set(handles.txtCalLevel,'string',70);
        set(handles.txtAnalysisStart,'Enable','off');
        set(handles.txtAnalysisEnd  ,'Enable','off');
        set(handles.txtOverlap,'enable','off');
        
        set(handles.chAvgLoudnessLimits,'enable','off');
        set(handles.chPercentiles,'enable','off');
        set(handles.chInternalNoise,'enable','on');
        set(handles.popInternalNoise,'enable','on');
        set(handles.popFc,'enable','on');
        set(handles.rbScripts,'value',1);
        set(handles.rbPsySound,'enable','off');
                
        set(handles.chParam1,'string','Low-passed filters');
        set(handles.chParam1,'enable','on');
        set(handles.chParam1,'value',1);
        
        set(handles.chParam2,'string','Auditory filterbank');
        set(handles.chParam2,'enable','on');
        set(handles.chParam2,'value',0);
        
        set(handles.chParam3,'string','Hair-cell envelope'); 
        set(handles.chParam3,'enable','off');
        
        set(handles.chParam4,'string','Adaptation loops');
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


% --- Executes on button press in chParam1.
function chParam1_Callback(hObject, eventdata, handles)
% hObject    handle to chParam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam1


% --- Executes on button press in chParam2.
function chParam2_Callback(hObject, eventdata, handles)
% hObject    handle to chParam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam2


% --- Executes on button press in chParam3.
function chParam3_Callback(hObject, eventdata, handles)
% hObject    handle to chParam3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam3


% --- Executes on button press in chParam4.
function chParam4_Callback(hObject, eventdata, handles)
% hObject    handle to chParam4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam4


% --- Executes on button press in chParam5.
function chParam5_Callback(hObject, eventdata, handles)
% hObject    handle to chParam5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam5


% --- Executes on button press in chParam6.
function chParam6_Callback(hObject, eventdata, handles)
% hObject    handle to chParam6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chParam6


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


% --- Executes on button press in chAvgLoudnessLimits.
function chAvgLoudnessLimits_Callback(hObject, eventdata, handles)
% hObject    handle to chAvgLoudnessLimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chAvgLoudnessLimits



function txtOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to txtOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOverlap as text
%        str2double(get(hObject,'String')) returns contents of txtOverlap as a double


% --- Executes during object creation, after setting all properties.
function txtOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popFc.
function popFc_Callback(hObject, eventdata, handles)
% hObject    handle to popFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popFc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popFc


% --- Executes during object creation, after setting all properties.
function popFc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popFc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles)
% hObject    handle to bSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSave



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


% --- Executes on button press in chFile1.
function chFile1_Callback(hObject, eventdata, handles)
% hObject    handle to chFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chFile1


% --- Executes on button press in chFile2.
function chFile2_Callback(hObject, eventdata, handles)
% hObject    handle to chFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chFile2



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X. Executes on button press in btnChange.
function btnChange_Callback(hObject, eventdata, handles)
%   - hObject    handle to btnChange (see GCBO)
%   - eventdata  reserved - to be defined in a future version of MATLAB
%   - handles    structure with handles and user data (see GUIDATA)

string1 = get(handles.txtFile1,'String');
string2 = get(handles.txtFile2,'String');
lbl1 = get(handles.txtLabel1,'String');
lbl2 = get(handles.txtLabel2,'String');

set(handles.txtFile1,'String',string2);
set(handles.txtFile2,'String',string1);
set(handles.txtLabel1,'String',lbl2);
set(handles.txtLabel2,'String',lbl1);

% --- Executes during object deletion, before destroying properties.
function btnCalculate_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to btnCalculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function popInternalNoise_Callback(hObject, eventdata, handles)
% hObject    handle to popInternalNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx = get(hObject,'Value');

if idx ~= 1
    set(handles.btnGetTemplate,'enable','on');
    set(handles.btnSimulateAFC,'enable','on');
else
    set(handles.btnGetTemplate,'enable','off');
    set(handles.btnSimulateAFC,'enable','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function popInternalNoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popInternalNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in chInternalNoise.
function chInternalNoise_Callback(hObject, eventdata, handles)
% hObject    handle to chInternalNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx = get(hObject,'Value'); % returns toggle state of chInternalNoise
if idx == 0
    set(handles.popInternalNoise,'Enable','off');
else
    set(handles.popInternalNoise,'Enable','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in chFile1.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to chFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chFile1


function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to txtGain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGain1 as text
%        str2double(get(hObject,'String')) returns contents of txtGain1 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popStepdB.
function popStepdB_Callback(hObject, eventdata, handles)
% hObject    handle to popStepdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popStepdB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popStepdB


% --- Executes during object creation, after setting all properties.
function popStepdB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popStepdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popNreversals.
function popNreversals_Callback(hObject, eventdata, handles)
% hObject    handle to popNreversals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popNreversals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popNreversals


% --- Executes during object creation, after setting all properties.
function popNreversals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popNreversals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chStepSizeHalved.
function chStepSizeHalved_Callback(hObject, eventdata, handles)
% hObject    handle to chStepSizeHalved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chStepSizeHalved


% --- Executes on button press in chAddDFT.
function chAddDFT_Callback(hObject, eventdata, handles)
% hObject    handle to chAddDFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chAddDFT

% --- Executes on selection change in popExamples.
function popExamples_Callback(hObject, eventdata, handles)
% hObject    handle to popExamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nExample = il_get_nAnalyser( handles.popExamples );

filename = AMTControl_Examples(nExample);

set(handles.txtFile1,'String',filename{1});
set(handles.txtFile2,'String',filename{2});
if nExample == 3
    handles.audio.filenameBBN = filename{3};
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popExamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popExamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popNavg.
function popNavg_Callback(hObject, eventdata, handles)
% hObject    handle to popNavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popNavg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popNavg


% --- Executes during object creation, after setting all properties.
function popNavg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popNavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popFc2.
function popFc2_Callback(hObject, eventdata, handles)
% hObject    handle to popFc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popFc2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popFc2


% --- Executes during object creation, after setting all properties.
function popFc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popFc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popStochasticDuration.
function popStochasticDuration_Callback(hObject, eventdata, handles)
% hObject    handle to popStochasticDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popStochasticDuration contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popStochasticDuration


% --- Executes during object creation, after setting all properties.
function popStochasticDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popStochasticDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popNsim.
function popNsim_Callback(hObject, eventdata, handles)
% hObject    handle to popNsim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popNsim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popNsim


% --- Executes during object creation, after setting all properties.
function popNsim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popNsim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline functions:
% 1 of 4:
function nAnalyser = il_get_nAnalyser(h)

strTmp = get(h,'String');
idxTmp = get(h,'value');
strTmp = strTmp{idxTmp};
strTmp = strsplit(strTmp,'-');
nAnalyser = str2num(strTmp{1});

function value = il_get_value_numericPop(h)

strTmp = get(h,'String');
idxTmp = get(h,'value');
strTmp = strTmp{idxTmp};
value = str2num(strTmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 of 4
function [bPlotParams, labelParams] = il_get_plotParams(handles,N)

if nargin < 2
    N = 7;
end

for i = 1:N
    exp1 = sprintf('bPlotParams(%.0f) = get(handles.chParam%.0f,''value''); labelParams{%.0f} = get(handles.chParam%.0f,''string'');',i,i,i,i);
    eval( exp1 );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 of 4:
function h = il_Plot_dau1996(stOuts,opts)

opts = ef(opts,'Title',''); 
opts = ef(opts,'bSave',0); % to be used in Mesh

out = stOuts.out;
fs  = stOuts.fs;
fc  = stOuts.fc; 
t = (1:size(out,1))/fs;

opts.step1  = 1;
opts.step2  = 1;
opts.bPlot3D= 0;
opts.bPlot2D= ~opts.bPlot3D;
opts.XLabel = 'Time [s]';
opts.YLabel = 'Centre frequency [Hz]';
opts.ZLabel = 'Amplitude [MU]';

if opts.bPlot2D
    opts.I_Matrix = [8,4];
    opts.I_Width  = opts.I_Matrix(2)*4;
    opts.I_Height = opts.I_Matrix(1)*4;
end

figure;
h = Mesh(t,fc,transpose(out),opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 of 4:
function h = il_Plot_dau1997(stOuts,opts)

opts = ef(opts,'Title',''); 
opts = ef(opts,'bSave',0); % to be used in Mesh

out = stOuts.out;
fs  = stOuts.fs;
fcm = stOuts.fcm;
t = (1:length(out))/fs;

NCh_mod = size(out,2);
fmod_num= 1:NCh_mod;
fmod    = fcm(fmod_num);

Zmin = min(min(out(:,fmod_num)));
Zmax = max(max(out(:,fmod_num)));

opts.step1  = 1;
opts.step2  = 1;
opts.bPlot3D= 0;
opts.bPlot2D= ~opts.bPlot3D;
opts.XLabel = 'Time [s]';
opts.YLabel = 'Modulation frequency [Hz]';
opts.ZLabel = 'Amplitude [MU]';

if opts.bPlot2D
    if NCh_mod <= 4
        opts.I_Matrix = [1,4];
    elseif NCh_mod <= 8
        opts.I_Matrix = [2,4];
    else
        opts.I_Matrix = [3,4];
    end
    opts.I_Width  = opts.I_Matrix(2)*4;
    opts.I_Height = opts.I_Matrix(1)*4;
end

figure;
h = Mesh(t,fmod,transpose(out(:,fmod_num)),opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_randomise_insig(insig)

Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];

function [outsig1 outsig2] = il_random_sample_mod_experiment(insig,fs,fc,BW,Fmod,ModIndex,dBFS)

if nargin < 7
    dBFS = 100;
end

N = length(insig);
r = cos_ramp(N,fs,200,200); r = r(:);
insig1 = il_randomise_insig(insig);
lvl = rmsdb(insig1) + dBFS;
insig1 = Set_Fourier_coeff_to_zero(insig1,fs,fc-BW/2,fc+BW/2);
insig1 = r .* insig1;
outsig1 = setdbspl(insig1,lvl,dBFS);

insig2 = il_randomise_insig(insig);
insig2 = Set_Fourier_coeff_to_zero(insig2,fs,fc-BW/2,fc+BW/2);
insig2 = r .* insig2;

[depthFastl depthDau] = m2d(ModIndex);
                
[insig2 env] = ch_am(insig2,Fmod,depthFastl,'d',fs,0);
outsig2 = setdbspl(insig2,lvl,dBFS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = il_pool_in_one_column(incell)

out = [];
for k = 1:length(incell);
    outtmp = incell{k};
    out = [out; outtmp(:)];
end
y = out;


% --- Executes on button press in chRampMasker.
function chRampMasker_Callback(hObject, eventdata, handles)
% hObject    handle to chRampMasker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chRampMasker


% --- Executes on button press in chRampSignal.
function chRampSignal_Callback(hObject, eventdata, handles)
% hObject    handle to chRampSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chRampSignal


% --- Executes on selection change in popGainSupra.
function popGainSupra_Callback(hObject, eventdata, handles)
% hObject    handle to popGainSupra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popGainSupra contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popGainSupra


% --- Executes during object creation, after setting all properties.
function popGainSupra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popGainSupra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chSaveTemplate.
function chSaveTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to chSaveTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chSaveTemplate
