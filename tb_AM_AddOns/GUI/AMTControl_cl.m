function outs = AMTControl_cl(handles_man)
% function outs = AMTControl_cl(handles_man)
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
    handles = il_SimulateAFC(handles);

    disp(['Threshold (rel. in dB): ' num2str(handles.Threshold)])
    disp(['Variance: ' num2str(handles.sigma)])
    try
        disp(['Script in template   : ' handles.script_template])
        disp(['Script in simulations: ' handles.script_sim])
    end

    outs.Threshold = handles.Threshold;
    outs.Staircase = handles.Staircase;
    outs.StaircaseCC = handles.StaircaseCC;
    outs.audio = handles.audio; % it includes template
end

%%% In September 2015, old results: 
%  jepsen2008preproc_1Ch
%         criterion
%         4    /   2
% -------------------
%   0.4   -1.5 /
%** 0.5   0.5  /
%   0.7   3.5  /
%   0.8   4.5  /
%   0.9   5.5  /

%  dau1996preproc_1Ch
%           criterion
%           4    /   2
% -------------------
%   0.5   -10.5  /
%   0.8    -6.5  /
%   1.1    -3.5  /
%   1.5    -0.5  /
%** 1.7     0.5  /


%%% On 24/09/2015: 
%  jepsen2008preproc_1Ch
%         criterion
%         4    /   2
% -------------------
%   0.4    /
%** 0.5     /
%   0.7     /
%   0.8     /
%   0.9     /

%  dau1996preproc_1Ch
%           criterion
%           4    /   2
% -------------------
%   0.5   0   /
%   0.8   2.5 /
%   1.1     /
%   1.5     /
%** 1.7     /


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
sigma   = 0;

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
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996a',sigma,Ntimes,fc2plot_idx,tmp);
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
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996',sigma,Ntimes,fc2plot_idx,tmp);
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
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1997',sigma,Ntimes,fc2plot_idx,tmp);
                
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
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-modfilterbank',sigma,Ntimes,fc2plot_idx,tmp);
                
        template = Get_template_append(out_1Mean,out_2Mean,fs_intrep); % figure; plot(out_2Mean-out_1Mean)
        handles.script_template = tmp.script_template;
        
        handles.script_template = 'jepsen2008preproc_multi'; 
        
	case 104
        
        tmp.resample_intrep = handles.resample_intrep;
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-lowpass',sigma,Ntimes,fc2plot_idx,tmp);
                
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
function handles = il_SimulateAFC(handles)

nAnalyser   = handles.nAnalyser;
Level_start = handles.Gain4supra;

Level_step_i    = handles.StepdB;
Reversals_stop  = handles.Nreversals;
bHalveStepSize  = handles.bStepSizeHalved;
Reversals4avg   = handles.Reversals4avg;
StepdBmin       = handles.StepdBmin;

fs = handles.audio.fs;
fs_intrep = handles.audio.fs_intrep;

Ntimes = handles.audio.Ntimes; %
Nsim   = handles.Nsim;
bListen = 0;

masker_interval_avg = [];

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

fc2plot_idx = handles.fc2plot_idx;
fc2plot_idx2 = handles.fc2plot_idx2;

%%%
if length(fc2plot_idx) == 1 & fc2plot_idx2 == fc2plot_idx;
    fc      = audtofreq(fc2plot_idx+2,'erb'); % used as input for single-channel modelling
elseif length(fc2plot_idx) == 1 & fc2plot_idx2 ~= fc2plot_idx;
    fcmin   = audtofreq(min(fc2plot_idx)+2,'erb'); 
    fcmax   = audtofreq(max(fc2plot_idx2)+2,'erb');
    fc2plot_idx = fc2plot_idx:fc2plot_idx2;
else
    fcmin   = audtofreq(min(fc2plot_idx)+2,'erb'); 
    fcmax   = audtofreq(max(fc2plot_idx)+2,'erb');
end

mu      = 0;   
bSigma = handles.bInternalNoise;
if bSigma
    sigma   = handles.sigma;
else
    sigma = 0;
end

if length(fc2plot_idx) == 1 & fc2plot_idx2 == fc2plot_idx;
    bSingleChannel = 1;
    bMultiChannel = 0;
    fc2plot_idx_1 = 1;
else
    bSingleChannel = 0;
    bMultiChannel = 1;
end
%%%

Threshold = [];
bUseRamp       = handles.bUseRamp;
bUseRampSignal = handles.bUseRampS;

if bUseRamp;       rampdn = handles.DurRamps; end % ramps in ms
if bUseRampSignal; rampsl = rampdn; end

switch handles.increment_method
    case 'level'
        if bUseRampSignal
            fprintf('Introducing %.0f-ms ramps into test signals\n',rampsl);
            insig2 = Do_cos_ramp(insig2,fs,rampsl);
        end
end

bAddSilence2noise = handles.bAddSilence2noise;

for k = 1:Nsim
    
    Nmaskers = 0;
    Level_step = Level_step_i;
    if k > 1
        fprintf('Running simulation: %.0f of %.0f. Last computed threshold=%.2f [dB, relative]\n',k,Nsim,Threshold(k-1));
    else
        fprintf('Running simulation: %.0f of %.0f.\n',k,Nsim);
    end
    
    Level_current   = Level_start;

    Staircase   = [];
    StaircaseCC = [];
    Reversals   = [];

    nWrong      = 0;
    nCorrect    = 0;
    nReversal   = -1;
    bSucceeded  = 1; % not failed
    
    while (nReversal < Reversals_stop) & (bSucceeded ==  1) % up to line 765

        N = length(insig2);
        if bDeterministic == 1
            insig1s0 = insig1;
            insig1s1 = insig1;
            insig1s2 = insig1;
        end
        if bStochastic == 1
            
            switch handles.increment_method
                case 'modulation-depth'
                    insig2 = il_randomise_insig(handles.audio.insig1orig);
                    insig2 = insig2(1:N); % we re-assign insig2
                    
                    if bUseRampSignal
                        if k == 1; fprintf('Introducing %.0f-ms ramps into test signals\n',rampsl); end;
                        insig2 = Do_cos_ramp(insig2,fs,rampsl);
                    end

            end
            
            insig1s0 = il_randomise_insig(handles.audio.insig1orig);
            insig1s1 = il_randomise_insig(handles.audio.insig1orig);
            insig1s2 = il_randomise_insig(handles.audio.insig1orig);
        end
        
        if bAddSilence2noise == 1
            Lno = round(handles.Silence2noise * fs);
            Nno = N - Lno;
        else
            Lno = 0;
            Nno = N;
        end
        
        insig1s0 = insig1s0( 1:Nno );
        insig1s1 = insig1s1( 1:Nno );
        insig1s2 = insig1s2( 1:Nno );
        
        if bUseRamp
            if k==1; fprintf('Introducing %.0f-ms ramps into maskers\n',rampdn); end
            insig1s0 = Do_cos_ramp(insig1s0, fs, rampdn);
            insig1s1 = Do_cos_ramp(insig1s1, fs, rampdn);
            insig1s2 = Do_cos_ramp(insig1s2, fs, rampdn);
        end
        
        insig1s0 = [insig1s0; Gen_silence(Lno/fs,fs)];
        insig1s1 = [insig1s1; Gen_silence(Lno/fs,fs)];
        insig1s2 = [insig1s2; Gen_silence(Lno/fs,fs)];
        
        Gain2apply  = From_dB(Level_current);
        
        switch handles.increment_method
            case 'level'
                insig2test  = Gain2apply * insig2;
                interval2 = insig1s1 + insig2test; % Noise + current signal
                
            case 'modulation-depth'
                mdepth = Gain2apply; % d2m(Gain2apply,'dau');
                fmod = handles.fmod;
                insig2test = Do_AM(insig2,fs,fmod,mdepth);
                interval2 = insig2test;
        end
        interval1   = insig1s0; % Only noise
        interval1s2 = insig1s2; % Only noise 2
        
        if bListen
            pause(2); sound(interval1, fs); pause(1); sound(interval1s2, fs);
            pause(1); sound(interval2, fs); pause(1);
        end
        
        tmp = [];
        tmp.masker_ramp_ms = 0; % ramps already applied
        tmp.bAddSilence2noise = 0; % already applied
        
        switch handles.MethodIntRep
            case 1
            
            if nAnalyser == 99;
                
                model = 'dau1996a';
                [out_interval1   xx fs_intrep] = Get_internalrep_stochastic(interval1  ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval1s2 xx fs_intrep] = Get_internalrep_stochastic(interval1s2,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval2   xx fs_intrep] = Get_internalrep_stochastic(interval2 ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                
                if bMultiChannel
                    handles.script_sim = 'dau1996apreproc';
                end
                
                if bSingleChannel
                    handles.script_sim = 'dau1996apreproc_1Ch';
                end
                
            elseif nAnalyser == 100;
                
                model = 'dau1996';
                
                [out_interval1   xx fs_intrep otmp] = Get_internalrep_stochastic(interval1  ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval1s2 xx fs_intrep] = Get_internalrep_stochastic(interval1s2,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval2   xx fs_intrep] = Get_internalrep_stochastic(interval2 ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                
                nchn_dec = length(fc2plot_idx);
                
                if bMultiChannel
                    handles.script_sim = 'dau1996preproc';
                end
                
                if bSingleChannel
                    handles.script_sim = 'dau1996preproc_1Ch';
                end

            elseif nAnalyser == 101

                model = 'dau1997';
                
                tmp.modfiltertype = handles.modfiltertype;
                tmp.resample_intrep = handles.resample_intrep;
                tmp.chn_modfilt = 1:12;%:12;
                
                [out_interval1   xx fs_intrep otmp] = Get_internalrep_stochastic(interval1,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval1s2 xx fs_intrep] = Get_internalrep_stochastic(interval1s2,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval2   xx fs_intrep] = Get_internalrep_stochastic(interval2 ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                
                mfc = otmp.mfc;
                nchn_dec = length(otmp.fc); % size(out_interval1,1) / (length(interval1) / (fs/fs_intrep));
                
                if bMultiChannel
                    handles.script_sim = 'dau1997preproc';
                end
                if bSingleChannel
                    handles.script_sim = 'dau1997preproc_1Ch';
                end

            elseif nAnalyser == 103
                
                model = 'jepsen2008'; % same as 'jepsen2008-modfilterbank'
                
                tmp.resample_intrep = handles.resample_intrep;
                [out_interval1   xx fs_intrep otmp] = Get_internalrep_stochastic(interval1,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval1s2 xx fs_intrep] = Get_internalrep_stochastic(interval1s2,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval2   xx fs_intrep] = Get_internalrep_stochastic(interval2 ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                
                nchn_dec = length(otmp.fc);
                
                if bMultiChannel
                    handles.script_sim = 'jepsen2008preproc_multi';
                end
                
                if bSingleChannel
                    handles.script_sim = 'jepsen2008preproc_1Ch';
                end

            elseif nAnalyser == 104
                
                model = 'jepsen2008-lowpass';
                
                tmp.resample_intrep = handles.resample_intrep;
                [out_interval1   xx fs_intrep otmp] = Get_internalrep_stochastic(interval1,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval1s2 xx fs_intrep] = Get_internalrep_stochastic(interval1s2,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                [out_interval2   xx fs_intrep] = Get_internalrep_stochastic(interval2 ,[],fs,model,0,Ntimes,fc2plot_idx,tmp);
                
                nchn_dec = length(otmp.fc);
                
                handles.script_sim = otmp.script_template; % either 'jepsen2008preproc_1Ch' or 'jepsen2008preproc_multi'
                
            end
            
        case 2
            if nAnalyser == 99.1
                modelc = 'dau1996apreproc';
            elseif nAnalyser == 100.1
                modelc = 'dau1996preproc';
            elseif nAnalyser == 101.1
                modelc = 'dau1997preproc';
            else
                error('')
            end
                
            out_interval1   = casprepresentation(interval1  ,modelc,{fs});
            out_interval1s2 = casprepresentation(interval1s2,modelc,{fs});
            out_interval2   = casprepresentation(interval2  ,modelc,{fs});
            
            out_interval1   = out_interval1(:,fc2plot_idx);
            out_interval1s2 = out_interval1s2(:,fc2plot_idx);
            out_interval2   = out_interval2(:,fc2plot_idx);
        end
        
        %%% Add internal noise:
        out_interval1   = Add_gaussian_noise(out_interval1,mu,sigma); % Add internal noise
        out_interval1s2 = Add_gaussian_noise(out_interval1s2,mu,sigma);
        out_interval2   = Add_gaussian_noise(out_interval2,mu,sigma); % Add internal noise

        [a b] = size(template);
        % one audio-frequency band but all the modulation filterbanks:

        sigint1     = reshape(out_interval1,a,b);
        sigint1s2   = reshape(out_interval1s2,a,b);
        sigint2     = reshape(out_interval2,a,b);
        %%%        

        intM = handles.audio.avg_intrep_M; % mean([out_N0 out_N1 out_N2],2);
        intrep_M0   = intM; 
        intrep_M1   = intM; 
        intrep_M2   = intM; 

        diff11 =   sigint1-intrep_M0;
        diff12 = sigint1s2-intrep_M1;
        diff20 =   sigint2-intrep_M2; %tt = 1:length(diff11); figure; plot(tt,diff11+30,tt,diff12,tt,diff20-30); legend('11','12','20')

        switch handles.bDecisionMethod 
            case{2,3,4}
                cc1 = optimaldetector(sigint2, sigint1);
                cc2 = optimaldetector(sigint2, sigint1s2);
                if cc2 > cc1
                    diffGreatestCC = sigint2-sigint1s2;
                else
                    diffGreatestCC = sigint2-sigint1;
                end
        end
        
        bDecisionMethod = handles.bDecisionMethod; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. First decision used:
        % The following was took out on 02/09/2015, because it is only
        % accounting the internal variability, so it is not completely correct:
        if bDecisionMethod == 1
            [decision(1) corrmue(:,1)] = optimaldetector(diff11,template);
            [decision(2) corrmue(:,2)] = optimaldetector(diff12,template);
            [decision(3) corrmue(:,3)] = optimaldetector(diff20,template);
            
            finaldecision = ( decision(3)-max(decision(1),decision(2)) )/sigma; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. Second decision used:
        if bDecisionMethod == 2
            switch handles.MethodIntRep
                case 1
                    
                    decision(1) = optimaldetector(diff11,template,fs_intrep);
                    decision(2) = optimaldetector(diff12,template,fs_intrep);
                    decision(3) = optimaldetector(diffGreatestCC,template,fs_intrep); 
                    
                    decision = decision / (nchn_dec); % number of audio channels
                    
                case 2
                    decision(1) = optimaldetector(diff11,template);
                    decision(2) = optimaldetector(diff12,template);
                    decision(3) = optimaldetector(diffGreatestCC,template);
            end
            
            % time4var = 10e-3;
            % samples4var = floor(time4var * fs_intrep);
            % varn = mean(std(buffer(diff20(:),samples4var))); % buffered in 2.5-ms sections
            
            [xx idx_decision] = max(decision);
            
            switch handles.MethodIntRep
                case 1
                    value_Tr = sigma; % varn; % sigma; % max(varn); 
                case 2
                    value_Tr = 1000; % arbitrary value
            end
                        
            value_De = decision(3);% -max(decision(1:2));
            
            if idx_decision == 3 & value_De > value_Tr
                idx_decision = 3;
            else
                idx_decision = 1; % randi(3);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. Fourth decision used:
        if bDecisionMethod == 3
            rule        = [1 2]; % 2-down, 1-up
            numint      = 3;
            Mtmp        = mean([diff11 diff12 diffGreatestCC]);
            Stmp        = sigma*[1 1 1];
            
            [bDecided(1) pr(1)] = ideal_observer( Mtmp(1), Stmp(1), numint, rule(2));
            [bDecided(2) pr(2)] = ideal_observer( Mtmp(2), Stmp(2), numint, rule(2));
            [bDecided(3) pr(3)] = ideal_observer( Mtmp(3), max(Stmp(1:2)), numint, rule(2));
            idx_tmp = find(bDecided == 1);
            [maxvalue, idx_decision] = max( pr );
            if bDecided( idx_decision )~= 1
                idx_decision = 1; % randi(3); % just random number
            end
            
            decision = pr(idx_decision);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. Fourth decision used:
        if bDecisionMethod == 4
            rule        = [1 2]; % 2-down, 1-up
            numint      = 3;
            % Mtmp        = abs(mean([diff11 diff12 diff20]));
            Mtmp        = abs( mean([diff11 diff12 diffGreatestCC]) );
            % Stmp        = std([diff11 diff12 diff20]);
            Stmp        = sigma*[1 1 1];
            [bDecided(1) pr(1)] = caspmdecide( Mtmp(1), Stmp(1), rule, numint);
            [bDecided(2) pr(2)] = caspmdecide( Mtmp(2), Stmp(2), rule, numint);
            [bDecided(3) pr(3)] = caspmdecide( Mtmp(3), max(Stmp(1:2)), rule, numint);
            % idx_tmp = find(bDecided == 1);
            [maxvalue, idx_decision] = max( pr );
            if bDecided( idx_decision )~= 1
                idx_decision = 1; % randi(3); % just random number
            end
            
            decision = pr(idx_decision);
        end
        Staircase = [Staircase; Level_current];
        StaircaseCC = [StaircaseCC decision'];
        
        if idx_decision == 3 
            if nCorrect >= 1

                if nWrong == 1
                    nReversal = nReversal + 1;
                    Reversals = [Reversals; Level_current];
                    if mod(nReversal,2) == 0 & bHalveStepSize & nReversal~= 0
                        Level_step = max( Level_step/2, StepdBmin );
                    end
                    nWrong = 0;
                end

                if mod(nCorrect,2) == 1 
                    Level_current = Level_current - Level_step; % we make it more difficult
                end
            else
                
                nWrong = 0;
            end 
            nCorrect = nCorrect+1;
            
        else % masker is chosen
            if nWrong == 0

                if nCorrect >= 1
                    nReversal = nReversal + 1;
                    Reversals = [Reversals; Level_current];
                    if mod(nReversal,2) == 0 & bHalveStepSize 
                        Level_step = max( Level_step/2, StepdBmin );
                    end
                end
                
                nCorrect = 0;
                
            end
            nWrong = nWrong+1;
            
            Level_current = Level_current + Level_step; % we make it easier after two mistakes
            if strcmp(handles.increment_method,'modulation-depth') & Level_current > 0
                warning('Modulation depth (dB) has to be less or equal than 0 dB...')
                Level_current = Level_current - Level_step;
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
        Threshold(k) = median(Reversals(end-Reversals4avg+1:end,:));
    else
        Threshold(k) = NaN;
    end
    
    if handles.bDebug
        if bDeterministic
            figure; 
            plot(Staircase,'o');
            xlabel('Presentation order')
            ylabel('Relative level')
            title(sprintf('Reversals median = %.2f [dB]',Threshold(k)));
            grid on
        end

        if bStochastic & k == Nsim
            if nargout == 0
                figure;
                plot(Threshold,'o');
                xlabel('Running noise nr.')
                ylabel('Level at threshold (relative level)')
                title(sprintf('Average threshold median (L25,L75) = %.2f dB (%.2f, %.2f)',prctile(Threshold,50),prctile(Threshold,25),prctile(Threshold,75)));
                grid on
            end
        end
    end
    
    masker_interval_avg = []; % make sure it is deleted
    
end
handles.Threshold = Threshold;
handles.Staircase = Staircase;
handles.StaircaseCC = StaircaseCC;

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
