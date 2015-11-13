function [outsig, fc, params] = drnl_CASP_debug(insig,fs,varargin)
% function [outsig, fc, params] = drnl_CASP_debug(insig,fs,varargin)
% 
% 1. Description:
%       Script for PEMO preprocessing using an implementation of the dual
%       resonance nonlinear (DRNL) filter (Lopez-Poveda, meddis 2001)
%       The filter models the BM non-linearity. The output outsig is assumed
%       to be in [m/s], representing the basilar membrane velocity.
% 
% Author: Morten Loeve Jepsen, 2.nov 2005
%
% usage: out = drnl(x,CF,fs)
% Original file: '..\Psychoacoustics\CASP_Jepsen\CreateIntRepV02\casp2008\bm\drnl.m'
% Created on    : 24/04/2015
% Last update on: 15/08/2015
% Last use on   : 27/10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables in script getDRNLparam.m were all replaced and validated:
%   linDRNLpar(1).vals (kv.lin_fc) 
%   linDRNLpar(2).vals (kv.lin_ngt)
%   linDRNLpar(3).vals
%   linDRNLpar(4).vals
%   linDRNLpar(5).vals (kv.lin_lp_cutoff)
%   linDRNLpar(6).vals

%   nlinDRNLpar(1).vals (kv.nlin_fc_before)
%   nlinDRNLpar(2).vals (kv.nlin_ngt_before). This value is 2 but in Lopez-Poveda is three
%                       (kv.nlin_ngt_after).
%   nlinDRNLpar(3).vals (kv.nlin_bw_before)
%   nlinDRNLpar(4).vals (kv.nlin_a)
%   nlinDRNLpar(5).vals (kv.nlin_b)
%   nlinDRNLpar(6).vals (kv.nlin_c)
%   nlinDRNLpar(7).vals (kv.nlin_lp_cutoff)

% if nargin < 3
%     % necessary parameters for basilar-membrane and modulation filterbank
%     [BM MF]     = CaspPreProcCfg;
%     [BM MF Lp]  = CaspPreProcInit(BM, MF, fs);
%     warning('In this implementation IntRep neither MF are not being used...')
% end

% Import the parameters from the arg_drnl.m function.
definput.import={'drnl_CASP'};

[flags,kv,flow,fhigh]=ltfatarghelper({'flow','fhigh'},definput,varargin);

% Obtain the dboffset currently used.
dboffset=dbspl(1); % Assuming convention of 100 dB = rms of 1

% Switch signal to the correct scaling.
insig=gaindb(insig,dboffset-94); % Scales insig to Pa (gain of 6 dB if convention of 100 dB is used)

%% 1. Headphone filter (outer ear)
%       Typical human headphone-to-eardrum gain. Input insig is assumed to 
%       be in [Pa]
if flags.do_outerear
    hp_fir = headphonefilter(fs);   % Getting the filter coefficients at fs
    Ntmp = ceil(length(hp_fir)/2);  % Added by AO
    insig = [insig; zeros(Ntmp,1)]; % Added by AO
    insig = filter(hp_fir,1,insig);
    insig = insig(Ntmp+1:end);      % Added by AO
end

%% 2. Middle-ear filter
if flags.do_jepsenmiddleear
    me_fir = middleearfilter(fs,'jepsenmiddleear');
    insig = filter(me_fir,1,insig);
end

if flags.do_middleear
    
    me_fir = middleearfilter(fs);
    insig = filter(me_fir,1,insig); % This is already a minimum phase filter 
    me_gain_TF = max( 20*log10(abs(freqz(me_fir,1,8192))) ); % equivalent to FFT of 8192*2 points
    params.me_gain_TF = me_gain_TF;
    
end

if  flags.do_nomiddleear % added by AO
    
    me_fir = middleearfilter(fs); % correction using middleear from Lopez-Poveda
    me_gain_TF = max( 20*log10(abs(freqz(me_fir,1,8192))) ); % equivalent to FFT of 8192*2 points
    insig = gaindb(insig,me_gain_TF);
    params.me_gain_TF = me_gain_TF;
    
end
params.insig_after_TF = insig;

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, kv.bwmul, kv.basef); % Exactly the same freqs as in BM.CenterFreqs;

siglen      = length(insig);
nchannels   = length(fc);
outsig      = zeros(siglen,nchannels);
outsiglin   = zeros(siglen,nchannels);
outsignlin  = zeros(siglen,nchannels);

params.lin_gain = [];
params.lin_bw   = [];

lin_ngt = kv.lin_ngt;
lin_nlp = kv.lin_nlp;

nlin_ngt_before= kv.nlin_ngt_before;
nlin_ngt_after = kv.nlin_ngt_after;
nlin_nlp       = kv.nlin_nlp;

for ii = 1:nchannels
    
    CF = fc(ii); % current centre frequency
    
    % Params to be used in the linear part:
    lin_gain        = polfun(     kv.lin_gain, fc(ii)); % Linear gain
    lin_fc          = polfun(       kv.lin_fc, fc(ii)); 
    lin_bw          = polfun(       kv.lin_bw, fc(ii));
    lin_lp_cutoff   = polfun(kv.lin_lp_cutoff, fc(ii));
    
    % Params to be used in the non-linear part:
    nlin_fc_before  = polfun(kv.nlin_fc_before, fc(ii));
    nlin_bw_before  = polfun(kv.nlin_bw_before,fc(ii));
    nlin_lp_cutoff  = polfun(kv.nlin_lp_cutoff, fc(ii));
    
    if CF <= 1500
        nlin_a      = polfun(kv.nlin_a,fc(ii));
        nlin_b      = polfun(kv.nlin_b,fc(ii));
        nlin_c      = polfun(kv.nlin_c,fc(ii));
    else
        fcabove = 1500;
        nlin_a      = polfun(kv.nlin_a_above,fcabove);
        nlin_b      = polfun(kv.nlin_b_above,fcabove);
        nlin_c      = polfun(kv.nlin_c,fcabove);
    end
    
    [GTlin_b,GTlin_a] = coefGtDRNL(lin_fc,lin_bw,fs); % get GT filter coeffs
    [LPlin_b,LPlin_a] = coefLPDRNL(lin_lp_cutoff,fs); % get LP filter coeffs
    
    params.lin_gain(end+1) = lin_gain;
    params.lin_bw(end+1)   = lin_bw;
    
    % -------------- linear part ------------------------------------------
    if flags.do_bothparts || flags.do_linonly
        
        y_lin = lin_gain.*insig; % Apply linear gain
        
        % Gammatone filtering multiple times (cascade):
        for n = 1:lin_ngt 
            y_lin = real(filter(GTlin_b,GTlin_a,y_lin));
        end
        
        % Low-pass filtering multiple times (cascade):
        for n = 1:lin_nlp % cascade of lowpass filters
            y_lin = filter(LPlin_b,LPlin_a,y_lin);
        end
        
    else
        y_lin = zeros(size(insig));
    end
    
    % -------------- Non-linear part ------------------------------
    if flags.do_bothparts || flags.do_nlinonly
    
        [GTnlin_b,GTnlin_a] = coefGtDRNL(nlin_fc_before,nlin_bw_before,fs); % get GT filter coeffs
        [LPnlin_b,LPnlin_a] = coefLPDRNL(nlin_lp_cutoff,fs); % get LP filter coeffs

        y_nlin = insig;

        % Gammatone filtering multiple times (cascade):
        for n = 1:nlin_ngt_before % Gammatone filtering multiple times for cascading
            y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));
        end

        % Broken stick nonlinearity
        y_decide = transpose( [nlin_a*abs(y_nlin) nlin_b*(abs(y_nlin)).^nlin_c] ); 

        y_nlin = transpose( transpose( sign(y_nlin) ).* min(y_decide) );

        % Now GT filtering again
        for n = 1:nlin_ngt_after % Gammatone filtering multiple times for cascading
            y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));
        end
        
        % then LP filtering
        for n = 1:nlin_nlp % cascade of lowpass filters
            y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);
        end
            
        params.fc_abc(ii,:)    = [CF nlin_a nlin_b nlin_c lin_gain];
    else
        y_nlin = zeros(size(insig));
    end
    
    if flags.do_output_gain
        GaindB=kv.gain_after_drnl; 
    end
    if flags.do_no_output_gain
        GaindB=0;
    end

    outsig(:,ii)        = From_dB( GaindB )* (y_lin + y_nlin);
    outsiglin(:,ii)     = From_dB( GaindB ) * y_lin;
    outsignlin(:,ii)    = From_dB( GaindB ) * y_nlin;
end

params.kv = kv;
params.outsiglin    = outsiglin;
params.outsignlin   = outsignlin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eof

% Inline functions:
function outpar=polfun(par,fc)
% p0 = par(1)
% m  = par(2)
% outpar=10^(par(1)+par(2)*log10(fc));
outpar=10^(par(1))*fc^par(2);
