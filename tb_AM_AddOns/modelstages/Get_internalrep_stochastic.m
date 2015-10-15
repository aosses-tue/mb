function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc,opts)
% function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc,opts)
% function [out_1avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,fs,model,sigma,Ntimes,idx_fc,opts)
%
% 1. Description:
%       Obtains an averaged internal representation of the signal in_signal_pre
%       superimposed on a random selected sample of the buffered noise in_masker_pre.
%       
%     Input parameters:
%       in_masker_pre - buffered 'noise'
%       in_signal_pre - current test signal (without noise)
% 
%     Output parameters:
%       out_1avg corresponds to the averaged internal representation (but
%           without internal noise) of the Masker alone interval.
%       out_2avg corresponds to the averaged internal representation (but
%           without internal noise) of the Masker plus signal interval.
% 
%           setup.fs - Sampling frequency: relevant for the normalisation process.
% 
% 2. Stand-alone example:
%           filename = AMTControl_Examples(0);
%           [insig1 fs]      = Wavread(filename{1}); 
%           [insig2supra fs] = Wavread(filename{2}); insig2supra = insig2supra(1:0.1*fs);
%           Get_internalrep_stochastic(insig1, insig2supra,fs,'dau1996',0.85,10);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/08/2015
% Last update on: 01/09/2015 
% Last use on   : 02/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
if nargin < 8
    opts = [];
end 

if nargin < 7
    fc = 1000;
    fctmp = ceil( freqtoaud(fc,'erb') );
    fc = audtofreq(fctmp,'erb');
else

    if length(idx_fc) == 1
        fc = audtofreq(idx_fc+3-1); % assuming default fc values
        idx_fc = 1; % output of the single channel will have only one column
    else
        fcmin = audtofreq(min(idx_fc)+2,'erb'); 
        fcmax = audtofreq(max(idx_fc)+2,'erb');
    end

end

if nargin < 6
    Ntimes = 10; % reassigned later in case no in_signal_pre is specified
end

if nargin < 5
    sigma = 0.8;
end
%%%

if length(idx_fc) == 1
    bSingleChannel = 1;
    bMultiChannel = 0;
else
    bSingleChannel = 0;
    bMultiChannel = 1;
end 

switch model
    case {'dau1997','jepsen2008','jepsen2008-modfilterbank'}
        opts = ef(opts,'chn_modfilt',1:12);
        opts = Ensure_field(opts,'resample_intrep', 'resample_intrep');
end

N = size(in_signal_pre,1);
if N ~= 0
    bDo_only_one_intrep = 0;
else
    N = size(in_masker_pre,1);
    bDo_only_one_intrep = 1;
end

if bDo_only_one_intrep == 1
    Ntimes = 1;
end
    
mu = 0;
bAvgMethod = 2; % more optimised average method (not generating column array)

if bAvgMethod == 2
    Nmaskers = 0;
end

out_1 = [];
out_2 = [];

opts = Ensure_field(opts,'masker_ramp_ms',0);
opts = ef(opts,'signal_ramp_ms',0);
opts = ef(opts,'increment_method','level'); % default
increment_method = opts.increment_method;

switch model
    case 'dau1997'
        opts = Ensure_field(opts,'modfiltertype','dau1997');
end

masker_ramp_ms = opts.masker_ramp_ms;
signal_ramp_ms = opts.signal_ramp_ms;

if Ntimes ~= 1
    bStochastic = 1;
else
    bStochastic = 0;
end
bDeterministic = ~bStochastic;

%%%
if bDo_only_one_intrep == 0
    if signal_ramp_ms ~= 0
        display('Introducing ramp for the input signal as well')
        in_signal = Do_cos_ramp(in_signal_pre,fs,signal_ramp_ms);
    else
        in_signal = in_signal_pre;
    end
end
%%%

bAddSilence2noise = opts.bAddSilence2noise;

for i = 1:Ntimes
    
    if bStochastic
        in_masker_s0 = Randomise_insig(in_masker_pre); % random sample of the noise
    end
    if bDeterministic
        in_masker_s0 = in_masker_pre; 
    end
    
    if bStochastic
        in_masker_s1 = Randomise_insig(in_masker_pre); % random sample of the noise
    end
    if bDeterministic
        in_masker_s1 = in_masker_pre; % random sample of the noise
    end
    
    if bDo_only_one_intrep == 0
        if bAddSilence2noise == 1
            Lno = round(opts.Silence2noise * fs);
            Nno = N - Lno;
        else
            Lno = 0;
            Nno = N;
        end
        in_masker_s0 = in_masker_s0(1:Nno);
        in_masker_s1 = in_masker_s1(1:Nno);
    else
        Lno = 0;
    end
    
    if masker_ramp_ms ~= 0 % introducing ramps into masker signals
        if i == 1; display('Introducing ramp for the maskers'); end
        in_masker_s0 = Do_cos_ramp( in_masker_s0,fs,masker_ramp_ms );
        in_masker_s1 = Do_cos_ramp( in_masker_s1,fs,masker_ramp_ms );
    end
    
    if bAddSilence2noise
        in_masker_s0 = [in_masker_s0; Gen_silence(Lno/fs,fs)];
        in_masker_s1 = [in_masker_s1; Gen_silence(Lno/fs,fs)];
    end
    
    intervalN0 = in_masker_s0;
    
    if bDo_only_one_intrep == 0
        switch increment_method
            case 'level' % default
                intervalSN = in_masker_s1 + in_signal;
            case 'modulation-depth'
                intervalSN = in_signal;
        end
    end
    
    switch model
        
        case 'dau1996a'
            
            if bMultiChannel == 1
                [out_pre , fc] = dau1996apreproc(intervalN0,fs);
                out_pre = out_pre(:,idx_fc);
            end
            if bSingleChannel == 1
                [out_pre , fc] = dau1996apreproc_1Ch(intervalN0,fs,fc); 
            end
            
            [Ni,Mi] = size(out_pre);
            fs_intrep = fs;
            
            if bAvgMethod == 1
            
                tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
                out_1 = [out_1 tmp(:)]; % out_1 affected by internal noise
                
            elseif bAvgMethod == 2
                
                Nmaskers_c = Nmaskers;
                Nmaskers = Nmaskers + 1;
                if Nmaskers_c > 0
                    out_1 = (out_1*Nmaskers_c + out_pre(:) )/Nmaskers;
                else
                    out_1 = (                   out_pre(:) )/Nmaskers;
                end
                
            end
            
            if bDo_only_one_intrep == 0
                if bMultiChannel == 1
                    [out_pre , fc] = dau1996apreproc(intervalSN,fs); % out_2pre affected by external noise
                    out_pre = out_pre(:,idx_fc);
                    outs.script_template = 'dau1996apreproc';
                end
                if bSingleChannel == 1
                    [out_pre , fc] = dau1996apreproc_1Ch(intervalSN,fs,fc); % out_1pre affected by external noise
                    outs.script_template = 'dau1996apreproc_1Ch';
                end

                if bAvgMethod == 1

                    tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
                    out_2 = [out_2 tmp(:)]; 

                elseif bAvgMethod == 2

                    if Nmaskers_c > 0
                        out_2 = (out_2*Nmaskers_c + out_pre(:) )/Nmaskers;
                    else
                        out_2 = (                   out_pre(:) )/Nmaskers;
                    end

                end
            end
            
        case 'dau1996'
            
            if bMultiChannel == 1
                [out_pre , fc] = dau1996preproc(intervalN0,fs);
                out_pre = out_pre(:,idx_fc);
            end
            if bSingleChannel == 1
                [out_pre , fc] = dau1996preproc_1Ch(intervalN0,fs,fc); 
            end
            
            [Ni,Mi] = size(out_pre);
            fs_intrep = fs;
            
            if bAvgMethod == 1
            
                tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
                out_1 = [out_1 tmp(:)]; % out_1 affected by internal noise
                
            elseif bAvgMethod == 2
                
                Nmaskers_c = Nmaskers;
                Nmaskers = Nmaskers + 1;
                if Nmaskers_c > 0
                    out_1 = (out_1*Nmaskers_c + out_pre(:) )/Nmaskers;
                else
                    out_1 = (                   out_pre(:) )/Nmaskers;
                end
                
            end
            
            if bDo_only_one_intrep == 0
                if bMultiChannel == 1
                    [out_pre , fc] = dau1996preproc(intervalSN,fs); % out_2pre affected by external noise
                    out_pre = out_pre(:,idx_fc);
                    outs.script_template = 'dau1996preproc';
                end
                if bSingleChannel == 1
                    [out_pre , fc] = dau1996preproc_1Ch(intervalSN,fs,fc); % out_1pre affected by external noise
                    outs.script_template = 'dau1996preproc_1Ch';
                end

                if bAvgMethod == 1

                    tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
                    out_2 = [out_2 tmp(:)]; 

                elseif bAvgMethod == 2

                    if Nmaskers_c > 0
                        out_2 = (out_2*Nmaskers_c + out_pre(:) )/Nmaskers;
                    else
                        out_2 = (                   out_pre(:) )/Nmaskers;
                    end

                end
            end
            
        case 'dau1997'
            
            if bMultiChannel == 1
                
                chn_modfilt = opts.chn_modfilt;
                tmp = []; 
                [out_1pre, fc, mfc, fs_intrep] = dau1997preproc_multi(intervalN0, fs,fcmin,fcmax,opts.modfiltertype, opts.resample_intrep);
                
                for j = 1:length(idx_fc) % default:  1:12
                    tmp = [tmp out_1pre{j}];
                end
                out_1pre = tmp;
                
                if bDo_only_one_intrep == 0
                    tmp = [];
                    [out_2pre , fc, mfc] = dau1997preproc_multi(intervalSN,fs,fcmin,fcmax,opts.modfiltertype, opts.resample_intrep);
                    for j = 1:length(idx_fc) % default:  1:12
                        tmp = [tmp out_2pre{j}];
                    end
                    out_2pre = tmp;
                end
                
                outs.script_template = 'dau1997preproc';
                    
            end
            
            if bSingleChannel == 1
                [out_1pre , fc, mfc, fs_intrep] = dau1997preproc_1Ch(intervalN0,fs,fc,opts.modfiltertype,opts.resample_intrep);
                
                if bDo_only_one_intrep == 0
                    [out_2pre , fc, mfc] = dau1997preproc_1Ch(intervalSN,fs,fc,opts.modfiltertype,opts.resample_intrep);
                end
                chn_modfilt = opts.chn_modfilt;
                if length(chn_modfilt) > length(mfc)
                    chn_modfilt = 1:length(mfc);
                end
                    
                out_1pre = [out_1pre(:,chn_modfilt)]; %out_1pre(:,2);out_1pre(:,2)];
                if bDo_only_one_intrep == 0
                    out_2pre = [out_2pre(:,chn_modfilt)]; %out_2pre(:,2);out_2pre(:,2)];
                end
                outs.script_template = 'dau1997preproc_1Ch';
            end
            
            [n m] = size(out_1pre);
            
            if bAvgMethod == 1
                
                out_1 = [out_1 Add_gaussian_noise(out_1pre(:),mu,sigma)]; 
                if bDo_only_one_intrep == 0
                    out_2 = [out_2 Add_gaussian_noise(out_2pre(:),mu,sigma)]; % Add internal noise
                end
                
            elseif bAvgMethod == 2
                
                Nmaskers_c = Nmaskers;
                Nmaskers = Nmaskers + 1;
                if Nmaskers_c > 0
                    out_1 = (out_1*Nmaskers_c + out_1pre(:) )/Nmaskers;
                    if bDo_only_one_intrep == 0
                        out_2 = (out_2*Nmaskers_c + out_2pre(:) )/Nmaskers;
                    end
                else
                    out_1 = (                   out_1pre(:) )/Nmaskers;
                    if bDo_only_one_intrep == 0
                        out_2 = (               out_2pre(:) )/Nmaskers;
                    end
                end
                
            end
            
        case {'jepsen2008','jepsen2008-modfilterbank','jepsen2008-lowpass'}
            
            switch model
                case {'jepsen2008','jepsen2008-modfilterbank'}
                    fbstyle = 'modfilterbank';
                case 'jepsen2008-lowpass'
                    fbstyle = 'lowpass';
            end
            
            if bMultiChannel == 1
                [out_1pre , fc, mfc, IntRep] = jepsen2008preproc_multi(intervalN0,fs,fcmin,fcmax,fbstyle,opts.resample_intrep);
                if bDo_only_one_intrep == 0
                    [out_2pre , fc] = jepsen2008preproc_multi(intervalSN,fs,fcmin,fcmax,fbstyle,opts.resample_intrep);
                end
                [Ni,Mi] = size(out_1pre{1});
                outs.script_template = 'jepsen2008preproc_multi';
            end
            if bSingleChannel == 1
                [out_1pre , fc, mfc, IntRep] = jepsen2008preproc_1Ch(intervalN0,fs,fc,fbstyle,opts.resample_intrep);
                if bDo_only_one_intrep == 0
                    [out_2pre , fc]          = jepsen2008preproc_1Ch(intervalSN,fs,fc,fbstyle,opts.resample_intrep);
                end
                outs.script_template = 'jepsen2008preproc_1Ch';
            end
            fs_intrep = IntRep.fs_intrep;
            
            out_1pre = il_pool_in_one_column(out_1pre);
            if bDo_only_one_intrep == 0
                out_2pre = il_pool_in_one_column(out_2pre);
            end
            Mi = length(idx_fc);
            Ni = round(length(out_1pre)/Mi);
            
            % 'Running' noise
            [n m] = size(out_1pre);
            
            if bAvgMethod == 1
                
                out_1 = [out_1 Add_gaussian_noise(out_1pre(:),mu,sigma)]; 
                if bDo_only_one_intrep == 0
                    out_2 = [out_2 Add_gaussian_noise(out_2pre(:),mu,sigma)]; % Add internal noise
                end
                
            elseif bAvgMethod == 2
                
                Nmaskers_c = Nmaskers;
                Nmaskers = Nmaskers + 1;
                if Nmaskers_c > 0
                    out_1 = (out_1*Nmaskers_c + out_1pre(:) )/Nmaskers;
                    if bDo_only_one_intrep == 0
                        out_2 = (out_2*Nmaskers_c + out_2pre(:) )/Nmaskers;
                    end
                else
                    out_1 = (                   out_1pre(:) )/Nmaskers;
                    if bDo_only_one_intrep == 0
                        out_2 = (                   out_2pre(:) )/Nmaskers;
                    end
                end
                
            end
            
        otherwise
            error('Model not added yet')
    end
    
end
   
if bAvgMethod == 1
    out_1avg = mean(out_1,2);
    if bDo_only_one_intrep == 0
        out_2avg = mean(out_2,2);
    else
        out_2avg = [];
    end
else
    out_1avg = out_1; % average should not be strictly necessary...
    if bDo_only_one_intrep == 0
        out_2avg = out_2;
    else
        out_2avg = [];
    end
end

try 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Reshape tested/non-tested
    %                   dau1996  dau1997  jepsen2008   jepsen2008-lp
    % bSingleChannel    
    % bMultiChannel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out_1avg = reshape(out_1avg,Ni,Mi);
    if bDo_only_one_intrep == 0
        out_2avg = reshape(out_2avg,Ni,Mi);
    end
end

outs.fc = fc;
switch model
    case 'dau1997'
        outs.mfc = mfc;
    case {'jepsen2008','jepsen2008-modfilterbank'}
        outs.mfc = mfc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = il_pool_in_one_column(incell)

out = [];
if iscell(incell)
    for k = 1:length(incell);
        outtmp = incell{k};
        out = [out; outtmp(:)];
    end
else
    out = incell(:);
end
y = out;
