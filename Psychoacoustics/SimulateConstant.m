function handles = SimulateConstant(handles)
% function handles = SimulateConstant(handles)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 07/11/2015
% Last update on: 07/11/2015 
% Last use on   : 07/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nAnalyser       = handles.nAnalyser;
handles = Ensure_field(handles,'Level_start',handles.Gain4supra);

Level_start     = handles.Level_start;

Level_step_i    = handles.StepdB;
Reversals_stop  = handles.Nreversals;
bHalveStepSize  = handles.bStepSizeHalved;
Reversals4avg   = handles.Reversals4avg;
StepdBmin       = handles.StepdBmin;
nDown           = handles.nDown;

fs              = handles.audio.fs;
fs_intrep       = handles.audio.fs_intrep;

Ntimes          = handles.audio.Ntimes; 
Nsim            = handles.Nsim;

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

bDetected = [];

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
        
    Level_current = Level_start;

    bSucceeded  = 1; % not failed
    
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
        error('Removed on 7/11/2015')
    end

    %%% Add internal noise:
    [out_interval1   xx] = Add_gaussian_noise(out_interval1,mu,sigma); % Add internal noise
    [out_interval1s2 yy] = Add_gaussian_noise(out_interval1s2,mu,sigma); % Add internal noise
    [out_interval2   zz] = Add_gaussian_noise(out_interval2,mu,sigma); % Add internal noise

    [a b] = size(template);
    % one audio-frequency band but all the modulation filterbanks:

    sigint1     = reshape(out_interval1,a,b);
    sigint1s2   = reshape(out_interval1s2,a,b);
    sigint2     = reshape(out_interval2,a,b);
    %%%        

    sizeT = size(template);
    noise1 = normrnd(mu,sigma,sizeT(1),sizeT(2));
    noise2 = normrnd(mu,sigma,sizeT(1),sizeT(2));
    noise3 = normrnd(mu,sigma,sizeT(1),sizeT(2));

    intM = handles.audio.avg_intrep_M; % mean([out_N0 out_N1 out_N2],2);
    intrep_M0   = intM + noise1; 
    intrep_M1   = intM + noise2; 
    intrep_M2   = intM + noise3; 

    diff11 =   sigint1-intrep_M0;
    diff12 = sigint1s2-intrep_M1;
    diff20 =   sigint2-intrep_M2; %tt = 1:length(diff11); figure; plot(tt,diff11+30,tt,diff12,tt,diff20-30); legend('11','12','20')

    bDecisionMethod = handles.bDecisionMethod; 

    switch bDecisionMethod 
        case 2
            cc1 = optimaldetector(sigint2, sigint1);
            cc2 = optimaldetector(sigint2, sigint1s2);
            if cc2 > cc1
                diffGreatestCC = sigint2-sigint1s2;
            else
                diffGreatestCC = sigint2-sigint1;
            end
        case {1, 3, 4}
            error('Removed on 7/11/2015')
        case 5
            diffGreatestCC = diff20;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Second decision used:
    %       This decision criterion is incorrect according to the discussion
    %       with AK on 4/11/2015

    if bDecisionMethod == 2
        switch handles.MethodIntRep
            case 1

                decision(1) = optimaldetector(diffGreatestCC,template,fs_intrep); 
                decision(2) = optimaldetector(diff11,template,fs_intrep);
                decision(3) = optimaldetector(diff12,template,fs_intrep);

                decision = decision / (nchn_dec); % number of audio channels

            case 2
                decision(1) = optimaldetector(diffGreatestCC,template);
                decision(2) = optimaldetector(diff11,template);
                decision(3) = optimaldetector(diff12,template);

        end

        decision = max(decision,normrnd(mu,sigma,1,3));

        switch nDown
            case 1
                [xx idx_decision] = max(decision(1:2)); % considering only one noise interval
            case 2
                [xx idx_decision] = max(decision);

        end

        switch handles.MethodIntRep
            case 1
                value_Tr = sigma; % varn; % sigma; % max(varn); 
        end

        value_De = decision(1); %max(decision(1:2)); 

        if idx_decision == 1 & value_De > value_Tr
            idx_decision = 1;
        else
            idx_decision = randi(3); % uniform distribution
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. Fifth decision used:
    if bDecisionMethod == 5
        switch handles.MethodIntRep
            case 1
                decision(1) = optimaldetector(diff20,template,fs_intrep); 
                decision(2) = optimaldetector(diff11,template,fs_intrep);
                decision(3) = optimaldetector(diff12,template,fs_intrep);

                decision = decision / (nchn_dec);
        end

        [xx idx_decision] = max(decision);

        decision = max(decision,normrnd(mu,sigma,1,3));

        switch nDown
            case 1
                [xx idx_decision] = max(decision(1:2));
            case 2
                [xx idx_decision] = max(decision);
        end


        if idx_decision ~= 1; disp(''); end % This line is just for debugging purposes...

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if idx_decision == 1 
        
        bDetected(end+1,1) = 1;
        
    else % masker is chosen
        
        bDetected(end+1,1) = 0;

    end 
   
    masker_interval_avg = []; % make sure it is deleted
    
end
handles.bDetected = bDetected;
handles.Percentage = 100*sum(bDetected)/length(bDetected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_randomise_insig(insig)

Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];
