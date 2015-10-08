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
% Last modified on  : 14/09/2015
% Last used on      : 14/09/2015 % Remember to check compatibility with template_PsySoundCL.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    handles_man = [];
end

handles_man = Ensure_field(handles_man,'Gain4supra',5);
handles_man = Ensure_field(handles_man,'nAnalyser',100); % 100 = dau1996
                                                         % 103 = jepsen2008
handles_man = Ensure_field(handles_man,'bDecisionMethod',2); % 2 = cross-correlation criterion
                                                             % 4 = m-AFC criterion
handles_man = Ensure_field(handles_man,'DurRamps',150); % ms

dir_out = Get_TUe_paths('outputs');

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

handles.increment_method = handles_man.increment_method;
if strcmp(handles.increment_method,'modulation-depth')
    handles.fmod    = handles_man.fmod;
    handles.dur_test= handles_man.dur_test;
end
handles.bDebug     = handles_man.bDebug;
handles.dir_out    = dir_out;
handles.audio      = handles_man.audio;
handles.Gain2file1 = handles_man.Gain2file1;
handles.Gain2file2 = handles_man.Gain2file2;

%%%
handles = il_btnLoad(filename1,filename2,handles);
%%%

handles.Gain4supra  = handles_man.Gain4supra;
handles.nAnalyser   = handles_man.nAnalyser;
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

dir_out = handles.dir_out;
Mkdir(dir_out); % creates folder in case it does not exist:

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

if bUseRamp;  tmp.masker_ramp_ms = handles.DurRamps; else; tmp.masker_ramp_ms = 0; end % ramp time in ms
if bUseRampS; tmp.signal_ramp_ms = handles.DurRamps; else; tmp.signal_ramp_ms = 0; end % ramp time in ms
        
if bDeterministic == 0 % then stochastic
    insig1 = handles.audio.insig1orig;
    
    L = length(insig2);
    insig2 = il_randomise_insig(insig1);
    insig2 = insig2(1:L); % we re-assign insig2
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
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
        handles.script_template = tmp.script_template;
        
    case 99.1
        
        [template_test2 out_2Mean out_1Mean] = casptemplate(insig2supra+insig1,insig1,'dau1996apreproc',{fs});
        template_test = template_test2(:,fc2plot_idx);
        
        fs_intrep = fs;
        handles.script_template = 'dau1996apreproc';
        
    case 100
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1996',sigma,Ntimes,fc2plot_idx,tmp);
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
        handles.script_template = tmp.script_template;
        
    case 100.1
        
        [template_test2 out_2Mean out_1Mean] = casptemplate(insig2supra+insig1,insig1,'dau1996preproc',{fs});
        template_test = template_test2(:,fc2plot_idx);
        
        fs_intrep = fs;
        handles.script_template = 'dau1996preproc';
        
    case 101 
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'dau1997',sigma,Ntimes,fc2plot_idx,tmp);
                
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        idxnan = find(isnan(out_1Mean)); out_1Mean(idxnan) = []; out_2Mean(idxnan) = [];
        
        handles.script_template = tmp.script_template;
        
    case 103
        
        [out_1Mean out_2Mean fs_intrep tmp] = Get_internalrep_stochastic(insig1,insig2supra,fs,'jepsen2008-modfilterbank',sigma,Ntimes,fc2plot_idx,tmp);
                
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        handles.script_template = tmp.script_template;
        
	case 104
        
        error('Continue here')
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
                handles.script_template = 'jepsen2008preproc_1Ch';
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
        
        template_test = Get_template_append(out_1Mean,out_2Mean,fs_intrep);
        
end

template = template_test;

t = ( 1:size(out_1Mean,1) )/fs_intrep;

switch nAnalyser
    case {99,99.1,100,100.1}
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))

    case 101
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))

    case {103,104}
        figure;
        plot(t,template); grid on
        xlabel(sprintf('Time [s]\nFollow the instructions in the command window to continue with the AFC simulation'))
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
    fc2plot_idx = 1;
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

        L = length(insig2);
        if bDeterministic == 1
            insig1s0 = insig1;
            insig1s1 = insig1;
            insig1s2 = insig1;
        end
        if bStochastic == 1
            
            switch handles.increment_method
                case 'modulation-depth'
                    insig2 = il_randomise_insig(handles.audio.insig1orig);
                    insig2 = insig2(1:L); % we re-assign insig2
                    
                    if bUseRampSignal
                        if k == 1; fprintf('Introducing %.0f-ms ramps into test signals\n',rampsl); end;
                        insig2 = Do_cos_ramp(insig2,fs,rampsl);
                    end

            end
            
            insig1s0 = il_randomise_insig(handles.audio.insig1orig);
            insig1s1 = il_randomise_insig(handles.audio.insig1orig);
            insig1s2 = il_randomise_insig(handles.audio.insig1orig);
        end
        
        insig1s0 = insig1s0( 1:L );
        insig1s1 = insig1s1( 1:L );
        insig1s2 = insig1s2( 1:L );
        
        if bUseRamp
            if k==1; fprintf('Introducing %.0f-ms ramps into maskers\n',rampdn); end
            insig1s0 = Do_cos_ramp(insig1s0, fs, rampdn);
            insig1s1 = Do_cos_ramp(insig1s1, fs, rampdn);
            insig1s2 = Do_cos_ramp(insig1s2, fs, rampdn);
        end
        
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
            pause(2)
            sound(interval1, fs)
            pause(1);
            sound(interval1s2, fs)
            pause(1)
            sound(interval2, fs)
            pause(1)
        end
             
        switch handles.MethodIntRep
            case 1
            
            if nAnalyser == 99;
                if bMultiChannel
                    
                    [out_interval1  , fc] = dau1996apreproc(interval1,fs);
                    [out_interval1s2, fc] = dau1996apreproc(interval1s2,fs);
                    [out_interval2  , fc] = dau1996apreproc(interval2,fs);
                     
                    out_interval1   = out_interval1(:,fc2plot_idx);
                    out_interval1s2 = out_interval1s2(:,fc2plot_idx);
                    out_interval2   = out_interval2(:,fc2plot_idx);
                     
                    handles.script_sim = 'dau1996apreproc';
                end
                if bSingleChannel
                    [out_interval1]   = dau1996apreproc_1Ch(interval1,fs,fc);
                    [out_interval1s2] = dau1996apreproc_1Ch(interval1s2,fs,fc);
                    [out_interval2]   = dau1996apreproc_1Ch(interval2,fs,fc);
                    
                    handles.script_sim = 'dau1996apreproc_1Ch';
                end
            elseif nAnalyser == 100;
                if bMultiChannel
                    [out_interval1  , fc] = dau1996preproc(interval1,fs);
                    [out_interval1s2, fc] = dau1996preproc(interval1s2,fs);
                    [out_interval2  , fc] = dau1996preproc(interval2,fs);
                    
                    out_interval1   = out_interval1(:,fc2plot_idx);
                    out_interval1s2 = out_interval1s2(:,fc2plot_idx);
                    out_interval2   = out_interval2(:,fc2plot_idx);

                    handles.script_sim = 'dau1996preproc';
                end
                if bSingleChannel
                    [out_interval1]   = dau1996preproc_1Ch(interval1,fs,fc);
                    [out_interval1s2] = dau1996preproc_1Ch(interval1s2,fs,fc);
                    [out_interval2]   = dau1996preproc_1Ch(interval2,fs,fc);
                    
                    handles.script_sim = 'dau1996preproc_1Ch';
                end

            elseif nAnalyser == 101

                if bMultiChannel
                    [out_interval1  ,fc] = dau1997preproc(interval1  ,fs);
                    [out_interval1s2,fc] = dau1997preproc(interval1s2,fs);
                    [out_interval2  ,fc] = dau1997preproc(interval2  ,fs);
                    
                    out_interval1   = il_pool_in_one_column(out_interval1(fc2plot_idx));
                    out_interval1s2 = il_pool_in_one_column(out_interval1s2(fc2plot_idx));
                    out_interval2   = il_pool_in_one_column(out_interval2(fc2plot_idx));
                                    
                    handles.script_sim = 'dau1997preproc';
                end
                if bSingleChannel
                    [out_interval1  ,fc, xx, outsfilt11] = dau1997preproc_1Ch(interval1  ,fs,fc);
                    [out_interval1s2,fc, xx, outsfilt12] = dau1997preproc_1Ch(interval1s2,fs,fc);
                    [out_interval2  ,fc, xx, outsfilt2] = dau1997preproc_1Ch(interval2  ,fs,fc);
                    
                    out_interval1   = out_interval1(:,:);
                    out_interval1s2 = out_interval1s2(:,:);
                    out_interval2   = out_interval2(:,:);
                    
                    handles.script_sim = 'dau1997preproc_1Ch';

                    bListenToFilter = 0;
                    if bListenToFilter
                        warning('Listen to filter...')
                        fprintf('Level current = %.0f dB, actual level noise = %.1f dB\n',Level_current,rmsdb(outsfilt11.out01_filterbank)+100)
                        sound(outsfilt11.out01_filterbank,fs)
                        pause(0.5)
                        sound(outsfilt12.out01_filterbank,fs)
                        pause(0.5)
                        fprintf('Level current = %.0f dB, actual level signal+noise = %.1f\n',Level_current,rmsdb(outsfilt2.out01_filterbank)+100)
                        sound(outsfilt2.out01_filterbank,fs)
                        disp('Press any key to continue...')
                        pause();
                    end
                end

            elseif nAnalyser == 103
                if bMultiChannel
                    [out_interval1] = jepsen2008preproc_multi(interval1,fs,fcmin,fcmax,'resample_intrep');
                    [out_interval1s2] = jepsen2008preproc_multi(interval1s2,fs,fcmin,fcmax,'resample_intrep');
                    [out_interval2] = jepsen2008preproc_multi(interval2,fs,fcmin,fcmax,'resample_intrep');

                    out_interval1   = il_pool_in_one_column(out_interval1);
                    out_interval1s2 = il_pool_in_one_column(out_interval1s2);
                    out_interval2   = il_pool_in_one_column(out_interval2);

                    handles.script_sim = 'jepsen2008preproc_multi';
                end
                if bSingleChannel
                    [out_interval1 xx xx outsfilt11]   = jepsen2008preproc_1Ch(interval1,fs,fc);
                    [out_interval1s2 xx xx outsfilt12] = jepsen2008preproc_1Ch(interval1s2,fs,fc);
                    [out_interval2 xx xx outsfilt2]    = jepsen2008preproc_1Ch(interval2,fs,fc);

                    handles.script_sim = 'jepsen2008preproc_1Ch';

                    bListenToFilter = 0;
                    if bListenToFilter
                        warning('Listen to filter')
                        fprintf('Level current = %.0f dB, actual level noise = %.1f dB\n',Level_current,rmsdb(outsfilt11.out_filterbank)+100)
                        sound(outsfilt11.out_filterbank,fs)
                        pause(0.5)
                        sound(outsfilt12.out_filterbank,fs)
                        pause(0.5)
                        fprintf('Level current = %.0f dB, actual level signal+noise = %.1f\n',Level_current,rmsdb(outsfilt2.out_filterbank)+100)
                        sound(outsfilt2.out_filterbank,fs)
                        disp('Press any key to continue...')
                        pause();

                    end

                    out_interval1   = il_pool_in_one_column(out_interval1);
                    out_interval1s2 = il_pool_in_one_column(out_interval1s2);
                    out_interval2   = il_pool_in_one_column(out_interval2);
                end

            elseif nAnalyser == 104
                if bMultiChannel
                    error('Not implemented yet (on 26/08/2015)');
                end
                if bSingleChannel
                    [out_interval1 xx xx outsfilt1] = jepsen2008preproc_1Ch(interval1,fs,fc,'lowpass');
                    [out_interval1s2] = jepsen2008preproc_1Ch(interval1s2,fs,fc,'lowpass');
                    [out_interval2 xx xx outsfilt2] = jepsen2008preproc_1Ch(interval2,fs,fc,'lowpass');

                    handles.script_sim = 'jepsen2008preproc_1Ch';
                end
            end
            
            [out_interval1 tmp] = Add_gaussian_noise(out_interval1,mu,sigma); % Add internal noise
            out_interval1s2 = Add_gaussian_noise(out_interval1s2,mu,sigma);
            out_interval2   = Add_gaussian_noise(out_interval2,mu,sigma); % Add internal noise
            
            [a b] = size(template);
            % one audio-frequency band but all the modulation filterbanks:

            sigint1     = reshape(out_interval1,a,b);
            sigint1s2   = reshape(out_interval1s2,a,b);
            sigint2     = reshape(out_interval2,a,b);
            %%%        

            intM = handles.audio.avg_intrep_M; % mean([out_N0 out_N1 out_N2],2);
            intrep_M0   = intM; % out_N0; % Add_gaussian_noise(out_interval1s2,mu,sigma);
            intrep_M1   = intM; % out_N1; % Add_gaussian_noise(out_interval1,mu,sigma);
            intrep_M2   = intM; % out_N2; % Add_gaussian_noise(out_interval1,mu,sigma);

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
                
            disp('');
            
        case 2
            
            out_interval1   = casprepresentation(interval1  ,'dau1996preproc',{fs});
            out_interval1s2 = casprepresentation(interval1s2,'dau1996preproc',{fs});
            out_interval2   = casprepresentation(interval2  ,'dau1996preproc',{fs});
        
            NN = 1;
            diff11 = [];
            diff12 = [];
            diff20 = [];

            for kk = 1:NN
                diff11(:,kk) = Add_gaussian_noise(out_interval1(:,fc2plot_idx2),mu,sigma); % Add internal noise
                diff12(:,kk) = Add_gaussian_noise(out_interval1s2(:,fc2plot_idx2),mu,sigma);
                diff20(:,kk) = Add_gaussian_noise(out_interval2(:,fc2plot_idx2),mu,sigma); % Add internal noise
            end
            diff11 = mean(diff11,2)-handles.xxir;
            diff12 = mean(diff12,2)-handles.xxir;
            diff20 = mean(diff20,2)-handles.xxir;
            
        end
        
        bDebug = 0;
        if bDebug == 1
            figure; 
            subplot(4,1,1); plot(interval2,'r'); hold on; plot(interval1); ha = gca;
            subplot(4,1,2); plot(out_interval2,'r'); hold on; plot(out_interval1); ha(end+1) = gca; 
            subplot(4,1,3); plot(diff20,'r'); hold on; plot(diff11); ha(end+1) = gca; 
            subplot(4,1,4); plot(out_interval2-out_interval1); % hold on; plot(diff11); ha(end+1) = gca; 
            title(sprintf('current level = %.1f',Level_current))
            disp('')
            close
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
                    % if fs == fs_intrep
                    decision(1) = optimaldetector(diff11,template,fs_intrep);
                    decision(2) = optimaldetector(diff12,template,fs_intrep);
                    decision(3) = optimaldetector(diffGreatestCC,template,fs_intrep); % diff20
                    % end
                case 2
                    error('continue')
                    [decision(1) corrmue(:,1)] = optimaldetector(diff11,template);
                    [decision(2) corrmue(:,2)] = optimaldetector(diff12,template);
                    [decision(3) corrmue(:,3)] = optimaldetector(diffGreatestCC,template);
            end
            diffs = [diff11 diff12 diff20];
     
            time4var = 10e-3;
            samples4var = floor(time4var * fs_intrep);
            
            varn = mean(std(buffer(diff20(:),samples4var))); % buffered in 2.5-ms sections
            % varn2 = mean(std(buffer(diff20(:),100)));
            
            [xx idx_decision] = max(decision);
            
            value_Tr = sigma; % varn; % sigma; % max(varn); 
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
