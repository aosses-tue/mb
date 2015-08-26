function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc)
% function [out_1avg out_2avg fs_intrep outs] = Get_internalrep_stochastic(in_masker_pre,in_signal_pre,fs,model,sigma,Ntimes,idx_fc)
%
% 1. Description:
%       Obtains an averaged internal representation of the signal in_signal_pre
%       superimposed on a random selected sample of the buffered noise in_masker_pre.
%       
%       in_masker_pre - buffered 'noise'
%       in_signal_pre - current test signal
% 
%           ir1 - related to 'Noise alone'
%           ir2 - related to 'Suprathreshold signal' (signal well above threshold)
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
% Last update on: 13/08/2015 
% Last use on   : 13/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    Ntimes = 10;
end

if nargin < 5
    sigma = 0.8;
end

if length(idx_fc) == 1
    bSingleChannel = 1;
    bMultiChannel = 0;
else
    bSingleChannel = 0;
    bMultiChannel = 1;
end 

N = size(in_signal_pre,1);
mu = 0;

out_1 = [];
out_2 = [];
rampdn = 20;
for i = 1:Ntimes
    
    in_masker_s0 = Randomise_insig(in_masker_pre); % random sample of the noise
    %in_masker_s0 = in_masker_s0(1:N,:);
    in_masker_s0 = Do_cos_ramp( in_masker_s0(1:N),fs,rampdn );
    in_masker_s1 = Randomise_insig(in_masker_pre); % random sample of the noise
    %in_masker_s1 = in_masker_s1(1:N,:);
    in_masker_s1 = Do_cos_ramp( in_masker_s1(1:N),fs,rampdn );
    switch model
        case 'dau1996'
            
            if bMultiChannel == 1
                [out_pre , fc] = dau1996preproc(in_masker_s0,fs);
            end
            if bSingleChannel == 1
                [out_pre , fc] = dau1996preproc_1Ch(in_masker_s0,fs,fc); 
            end
            
            [Ni,Mi] = size(out_pre(:,idx_fc));
            fs_intrep = fs;
            
            tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
            out_1 = [out_1 tmp(:)]; % out_1 affected by internal noise
            
            if bMultiChannel == 1
                [out_pre , fc] = dau1996preproc(in_masker_s1 + in_signal_pre,fs); % out_2pre affected by external noise
            end
            if bSingleChannel == 1
                [out_pre , fc] = dau1996preproc_1Ch(in_masker_s1 + in_signal_pre,fs,fc); % out_1pre affected by external noise
            end
            
            tmp = Add_gaussian_noise(out_pre(:,idx_fc),mu,sigma);
            out_2 = [out_2 tmp(:)]; % out_1 affected by internal noise
            
        case 'dau1997'
            
            if bMultiChannel == 1
                error('Not implemented yet, continue here');
            end
            if bSingleChannel == 1
                [out_1pre , fc, mfc] = dau1997preproc_1Ch(in_masker_s0                ,fs,fc);
                [out_2pre , fc, mfc] = dau1997preproc_1Ch(in_masker_s1 + in_signal_pre,fs,fc);
                fs_intrep = fs;
            end
            
            % 'Running' noise
            [n m] = size(out_1pre);
            out_1 = [out_1 Add_gaussian_noise(out_1pre(:),mu,sigma)]; 
            out_2 = [out_2 Add_gaussian_noise(out_2pre(:),mu,sigma)]; % Add internal noise
            
        case 'jepsen2008'
            
            fbstyle = 'modfilterbank';
            if bMultiChannel == 1
                [out_1pre , fc, xx, IntRep] = jepsen2008preproc_multi(in_masker_s0  ,fs,fcmin,fcmax,fbstyle,'resample_intrep');
                [out_2pre , fc] = jepsen2008preproc_multi(in_masker_s1 + in_signal_pre,fs,fcmin,fcmax,fbstyle,'resample_intrep');
                [Ni,Mi] = size(out_1pre{1});
            end
            if bSingleChannel == 1
                [out_1pre , fc, xx, IntRep] = jepsen2008preproc_1Ch(in_masker_s0  ,fs,fc,fbstyle);
                [out_2pre , fc] = jepsen2008preproc_1Ch(in_masker_s1 + in_signal_pre,fs,fc,fbstyle);
            end
            fs_intrep = IntRep.fs_intrep;
            
            out_1pre = il_pool_in_one_column(out_1pre);
            out_2pre = il_pool_in_one_column(out_2pre);
            Mi = length(idx_fc);
            Ni = round(length(out_1pre)/Mi);
            
            % 'Running' noise
            [n m] = size(out_1pre);
            out_1 = [out_1 Add_gaussian_noise(out_1pre(:),mu,sigma)]; 
            out_2 = [out_2 Add_gaussian_noise(out_2pre(:),mu,sigma)]; % Add internal noise
            
        otherwise
            error('Model not added yet')
    end
    
end
   
out_1avg = mean(out_1,2);
out_2avg = mean(out_2,2);

try
    out_1avg = reshape(out_1avg,Ni,Mi);
    out_2avg = reshape(out_2avg,Ni,Mi);
catch
    warning('Template not reshaped, probably you are using the modulation filterbank and then not every filter has the same amount of elements');
end

outs.fc = fc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = il_pool_in_one_column(incell)

out = [];
for k = 1:length(incell);
    outtmp = incell{k};
    out = [out; outtmp(:)];
end
y = out;
