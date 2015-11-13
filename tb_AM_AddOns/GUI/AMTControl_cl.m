function [outs outsIR] = AMTControl_cl(handles_man)
% function [outs outsIR] = AMTControl_cl(handles_man)
%
% 1. Description: 
%       AMTCONTROL MATLAB code for AMTControl.fig
%      
%       AMTCONTROL, by itself, creates a new AMTCONTROL or raises the 
%       existing singleton*.
%
%       Line    Stage                               Last updated on
%       101     3. btnLoad                          07/08/2015
%       338     4. btnGetTemplate                   07/08/2015
%       547     5. btnSimulateAFC                   02/09/2015
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
% Created on        : 14/09/2015 (Snapshot of AMTControl.m)
% Last modified on  : 30/10/2015
% Last used on      : 30/10/2015 % Remember to check compatibility with template_PsySoundCL.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    handles_man = [];
end

handles_man = Ensure_field(handles_man,'Gain4supra',5);
handles_man = Ensure_field(handles_man,'nAnalyser',100); % 100 = dau1996
                                                         % 103 = jepsen2008
handles_man = Ensure_field(handles_man,'bDecisionMethod',2); % 2 = cross-correlation criterion
                                                             % 4 = m-AFC criterion

% some additional settings for audio files:
handles_man = Ensure_field(handles_man,'DurRamps',150); % ms
handles_man = Ensure_field(handles_man,'bAddSilence2noise',0);
if handles_man.bAddSilence2noise == 1
    handles_man = Ensure_field(handles_man,'Silence2noise',0); % s
end

handles_man = Ensure_field(handles_man,'increment_method','level'); % 'modulation-depth'
handles_man = Ensure_field(handles_man,'filename1',[Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-60-dB-dur-800-ms.wav']);
switch handles_man.increment_method
    case 'level'
        handles_man = Ensure_field(handles_man,'filename2',[Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-42-dB-dur-800-ms.wav']);
    case 'modulation-depth'
        handles_man = Ensure_field(handles_man,'filename2',handles_man.filename1);
        handles_man = Ensure_field(handles_man,'fmod',5); % Hz
        handles_man = Ensure_field(handles_man,'dur_test',1); % s
end
handles_man = Ensure_field(handles_man,'do_template'  , 1);
handles_man = Ensure_field(handles_man,'do_simulation', 1);
handles_man = Ensure_field(handles_man,'MethodIntRep' , 1); % my method to obtain internal representations
handles_man = Ensure_field(handles_man,'Nreversals'   ,12);
handles_man = Ensure_field(handles_man,'nDown'        , 2); % 2-down, 1-up
handles_man = Ensure_field(handles_man,'Ntimes'       , 1); % deterministic

handles_man = Ensure_field(handles_man,'StepdB',2); 
handles_man = Ensure_field(handles_man,'StepdBmin',1); 
handles_man = Ensure_field(handles_man,'Reversals4avg',6); 
handles_man = Ensure_field(handles_man,'Gain2file1',0);
handles_man = Ensure_field(handles_man,'Gain2file2',0);
handles_man = Ensure_field(handles_man,'fc2plot_idx' ,ceil( freqtoaud(1000,'erb') )-2);
handles_man = Ensure_field(handles_man,'fc2plot_idx2',ceil( freqtoaud(1000,'erb') )-2);
handles_man = Ensure_field(handles_man,'Nsim',1);
filename1 = handles_man.filename1;
filename2 = handles_man.filename2;

switch handles_man.nAnalyser
    case 101
        handles_man = Ensure_field(handles_man,'modfiltertype','dau1997'); % dau1997wLP, derleth2000, jepsen2008
        handles_man = Ensure_field(handles_man,'resample_intrep','resample_intrep');
    case {103,104}
        handles_man = Ensure_field(handles_man,'resample_intrep','resample_intrep');
end

handles = handles_man; % we pass all the parameters in handles_man to handles

%%%
handles = il_btnLoad(filename1,filename2,handles);
%%%

handles.Gain4supra  = handles_man.Gain4supra;
handles.nAnalyser   = handles_man.nAnalyser;
switch handles.nAnalyser
    case 101
        handles.modfiltertype = handles_man.modfiltertype;
        handles.resample_intrep = handles_man.resample_intrep;
    case 103
        handles.resample_intrep = handles_man.resample_intrep;
end

handles.fc2plot_idx = handles_man.fc2plot_idx;
handles.fc2plot_idx2= handles_man.fc2plot_idx2;
handles.Ntimes      = handles_man.Ntimes; % 1 == deterministic
handles.bUseRamp    = handles_man.bUseRamp;
handles.bUseRampS   = handles_man.bUseRampS;
handles.DurRamps    = handles_man.DurRamps; % ms

%%%
if handles_man.do_template == 1
    handles = il_GetTemplate(handles);
else
    handles.audio.template = handles_man.audio.template;
end
%%%

if handles_man.do_simulation == 1
    handles.StepdB          = handles_man.StepdB;
    handles.StepdBmin       = handles_man.StepdBmin;
    handles.Nreversals      = handles_man.Nreversals;
    handles.bStepSizeHalved = 1;
    handles.Reversals4avg   = handles_man.Reversals4avg;
    handles.Nsim            = handles_man.Nsim;
    handles.bInternalNoise  = 1;
    handles.sigma             = handles_man.sigma;
    handles.bDecisionMethod = handles_man.bDecisionMethod;
    handles.MethodIntRep    = handles_man.MethodIntRep;
    handles = il_SimulateExperiment(handles);

    try
        disp(['Threshold (rel. in dB): ' num2str(handles.Threshold)])
        disp(['Variance: ' num2str(handles.sigma)])
    end
    
    try
        disp(['Script in template   : ' handles.script_template])
        disp(['Script in simulations: ' handles.script_sim])
    end

    try
        outs.Threshold = handles.Threshold;
        outs.Staircase = handles.Staircase;
        outs.StaircaseCC = handles.StaircaseCC;
    catch
        outs.bDetected = handles.bDetected;
        outs.Percentage = handles.Percentage;
    end
    outs.audio = handles.audio; % it includes template
    
    if nargout > 1
        outsIR.out_SN_interval = handles.out_SN_interval;
        outsIR.out_N_interval = handles.out_N_interval;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Executes on button press in btnLoad.
function handles = il_btnLoad(filename1, filename2, handles)

G1 = handles.Gain2file1; % gain in dB for audio file 1
G2 = handles.Gain2file2; % gain in dB for audio file 2

[insig1,fs]  = Wavread(filename1);
[insig2,fs2] = Wavread(filename2);

if fs ~= fs2
    error('Two audio files should have the same sampling frequency');
end

if strcmp(handles.increment_method,'modulation-depth')
    try
        insig2 = insig2(1:round(handles.dur_test*fs));
    catch
        error('dur_test field is longer than the length of filename1, please re-assign this parameter and re-run the code');
    end
end

handles.audio = Ensure_field(handles.audio,'fs',fs);

if fs ~= handles.audio.fs
    
    insig1 = resample(insig1,handles.audio.fs,fs);
    insig2 = resample(insig2,handles.audio.fs,fs);
    fs = handles.audio.fs;
end

t1 = ( 0:length(insig1)-1 )/fs;
t2 = ( 0:length(insig2)-1 )/fs;

ti_samples  = 1; % handles.audio.ti_samples;
tf_samples  = min(length(t1),length(t2)); % handles.audio.tf_samples;

handles.audio.ti_samples = ti_samples;
handles.audio.tf_samples = tf_samples;

xliminf = ti_samples/fs;
xlimsup = tf_samples/fs;

if handles.bDebug
    figure;
    subplot(2,1,1)
    plot(t1,From_dB(G1)*insig1);
    xlim([xliminf xlimsup])

    subplot(2,1,2)
    plot(t2,From_dB(G2)*insig2,'r');
    xlim([xliminf xlimsup])
end

adjustmentvalue = 0; % values calibrated to 100 dB RMS = 0 dBFS

insig1_orig = From_dB(adjustmentvalue+G1) * insig1;    
insig2_orig = From_dB(adjustmentvalue+G2) * insig2;

insig1 = insig1_orig(ti_samples:tf_samples);
insig2 = insig2_orig(ti_samples:tf_samples);

handles.audio.filename1 = filename1;
handles.audio.filename2 = filename2;

dBFS = 100;
RMS1 = rmsdb(insig1) + dBFS; % Zwicker's correction
RMS2 = rmsdb(insig2) + dBFS; 

thres_silence = 1/3; % one-third of the median value of the envelope 
[xx xx xx RMS1nosil] = Rmssilence(insig1,fs,thres_silence);
[xx xx xx RMS2nosil] = Rmssilence(insig2,fs,thres_silence);

fprintf('RMS1 = %.1f [dB SPL] (%.1f no silent)\n',RMS1,RMS1nosil + dBFS);
fprintf('RMS2 = %.1f [dB SPL] (%.1f no silent)\n',RMS2,RMS2nosil + dBFS);

bDFT = handles.bDebug;

if bDFT
    windowtype = 'hanning';
    
    K  = length(insig1)/2;
    [xx y1dB f] = freqfft2(insig1,K,fs,windowtype,dBFS);   
    [xx y2dB  ] = freqfft2(insig2,K,fs,windowtype,dBFS);
    
    figure;
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

handles.audio.insig1orig = insig1_orig;     
handles.audio.insig1 = insig1;
handles.audio.insig2 = insig2;
handles.audio.G1 = G1;
handles.audio.G2 = G2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Executes on button press in btnGetTemplate.
function handles = il_GetTemplate(handles)

insig2 = handles.audio.insig2;
fs = handles.audio.fs;

Gain4supra = handles.Gain4supra;
increment_method = handles.increment_method;
nAnalyser  = handles.nAnalyser;
   
fc2plot_idx = handles.fc2plot_idx;
fc2plot_idx2= handles.fc2plot_idx2;

fc_ERB      = fc2plot_idx+2;

fc2plot_idx = fc2plot_idx:fc2plot_idx2;

if length(fc2plot_idx) == 1
    fc      = audtofreq(fc_ERB,'erb'); % used as input for single-channel modelling
else
    fcmin   = audtofreq(min(fc2plot_idx)+2,'erb'); 
    fcmax   = audtofreq(max(fc2plot_idx)+2,'erb');
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
    
template_test = [];
Ntimes = handles.Ntimes;
if Ntimes == 1
    bDeterministic = 1;
else
    bDeterministic = 0;
end

bUseRamp  = handles.bUseRamp;
bUseRampS = handles.bUseRampS;

out_1 = [];
out_2 = [];

mu    = 0;
sigmaT = 0; % sigma template

tmp.bAddSilence2noise   = handles.bAddSilence2noise;
if tmp.bAddSilence2noise
    tmp.Silence2noise   = handles.Silence2noise;
end

if bUseRamp;  tmp.masker_ramp_ms = handles.DurRamps; else; tmp.masker_ramp_ms = 0; end % ramp time in ms
if bUseRampS; tmp.signal_ramp_ms = handles.DurRamps; else; tmp.signal_ramp_ms = 0; end % ramp time in ms
        
if bDeterministic == 0 % then stochastic
    insig1 = handles.audio.insig1orig;
else
    insig1 = handles.audio.insig1;
end

switch increment_method
    case 'level' % default case
        insig2supra = From_dB(Gain4supra) * insig2;
    case 'modulation-depth'
        mdepth = d2m(Gain4supra,'dau');
        fmod = handles.fmod;
                
        insig2supra = Do_AM(insig2,fs,fmod,mdepth);
        tmp.increment_method = 'modulation-depth';
end
    

switch nAnalyser
    
    case 99
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996a',sigmaT,Ntimes,fc2plot_idx,tmp);
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
        handles.script_template = tmp.script_template;
        
    case 99.1
        
        [template_test2 out_2Mean out_1Mean] = casptemplate(insig2supra+insig1,insig1,'dau1996apreproc',{fs});
        template = template_test2(:,fc2plot_idx);
        template=template/rms(template(:));
        out_1Mean = out_1Mean(:,fc2plot_idx);
        out_2Mean = out_2Mean(:,fc2plot_idx);
        
        fs_intrep = fs;
        handles.script_template = 'dau1996apreproc';
        
    case 100
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996',sigmaT,Ntimes,fc2plot_idx,tmp);
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
        handles.script_template = tmp.script_template;
        
    case 100.1
        
        [template_test2 out_2Mean out_1Mean] = casptemplate(insig2supra+insig1,insig1,'dau1996preproc',{fs});
        template = template_test2(:,fc2plot_idx);
        template=template/rms(template(:));
        out_1Mean = out_1Mean(:,fc2plot_idx);
        out_2Mean = out_2Mean(:,fc2plot_idx);
        
        fs_intrep = fs;
        handles.script_template = 'dau1996preproc';
        
    case 101 
        
        tmp.modfiltertype   = handles.modfiltertype;
        tmp.chn_modfilt     = 1:12;%:12;
        tmp.resample_intrep = handles.resample_intrep;
        handles.chn_modfilt = tmp.chn_modfilt;
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1997',sigmaT,Ntimes,fc2plot_idx,tmp);
                
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        idxnan = find(isnan(out_1Mean)); out_1Mean(idxnan) = 0; out_2Mean(idxnan) = 0;
        
        handles.script_template = tmp.script_template;
    
    case 101.1
        
        [template_test2 out_2Mean out_1Mean] = casptemplate(insig2supra+insig1,insig1,'dau1997preproc',{fs});
        template = template_test2(:,fc2plot_idx);
        template=template/rms(template(:));
        out_1Mean = out1_Mean(:,fc2plot_idx);
        out_2Mean = out2_Mean(:,fc2plot_idx);
        
        fs_intrep = fs;
        handles.script_template = 'dau1997preproc';    
        
    case 103
        
        tmp.resample_intrep = handles.resample_intrep;
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-modfilterbank',sigmaT,Ntimes,fc2plot_idx,tmp);
                
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep); % figure; plot(out_2Mean-out_1Mean)
        handles.script_template = tmp.script_template;
        
        handles.script_template = 'jepsen2008preproc_multi'; 
        
	case 104
        
        tmp.resample_intrep = handles.resample_intrep;
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-lowpass',sigmaT,Ntimes,fc2plot_idx,tmp);
                
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep); % figure; plot(out_2Mean-out_1Mean)
        handles.script_template = tmp.script_template;
        
        handles.script_template = 'jepsen2008preproc_multi'; 
        
end

t = ( 1:size(template,1) )/fs_intrep;

%%%
if handles.bDebug
    figure;
    subplot(2,1,1)
    plot(t,out_1Mean); grid on; hold on
    plot(t,out_2Mean,'r');
    legend('M','MTc')
    ha = gca; 
    
    subplot(2,1,2)
    plot(t,out_2Mean-out_1Mean); grid on;
    ha(end+1) = gca;
    linkaxes(ha,'x');
    
    p = Get_date;
    Saveas(gcf,['fig-template-' p.date2print],'epsc');
    close
end
%%%

try
switch nAnalyser
    case {99,99.1,100,100.1}
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))

    case 101
        figure;
        plot(t,template(:)); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))

    case {103,104}
        figure;
        plot(t,template(:)); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))
end
end
handles.audio.template      = template;
handles.audio.bDeterministic= bDeterministic;
handles.audio.avg_intrep_M  = out_1Mean;
handles.audio.avg_intrep_MT = out_2Mean;
handles.audio.Ntimes        = Ntimes;
handles.audio.fs_intrep     = fs_intrep;

bListen         = 0;

if bListen
    sound(tmp.inM, fs)
    pause(2);
    sound(tmp.inMT, fs)
    pause(2)
    sound(tmp.inM, fs)
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Executes on button press in btnSimulateAFC.
function handles = il_SimulateExperiment(handles)

handles = Ensure_field(handles,'experiment_type','AFC');
experiment_type = handles.experiment_type;

switch experiment_type
    case 'AFC'
        handles = SimulateAFC(handles);
    case 'constant'
        handles = SimulateConstant(handles);
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
