function [outsig, fc, outs] = kohlrausch1992preproc(insig, fs, varargin);
% function [outsig, fc, outs] = kohlrausch1992preproc(insig, fs, varargin);
%
% 1. Description:
%       Auditory model from Kohlrausch et. al. 1992
% 
%   Usage: [outsig, fc, outs] = kohlrausch1992preproc(insig,fs);
%          [outsig, fc, outs] = kohlrausch1992preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%     extra  : struct with outputs of each stage of the model (added by AO)
%
%   *Warning:* This code is incorrect, the Kohlrausch 1992 models uses a transmission-line 
%   model from Strube 1985, while THIS code erroneously uses the gammatone filters. 
%   If/when the Strube model is included in AMToolbox, this function will be fixed.
%
%   KOHLRAUSCH1992PREPROC(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz as described in Kohlrausch,
%   Puschel and Alphei (1992).
%  
%   [outsig,fc]=KOHLRAUSCH1992PREPROC(...) additionally returns the centre 
%   frequencies of the filter bank.
%
%   The Kohlrausch 1992 model consists of the following stages:
%   
%     1) a gammatone filter bank with 1-ERB spaced filtes.
%
%     2) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     3) an adaptation stage modelling nerve adaptation by a cascade of 5 loops.
%
%     4) a modulation low pass filter liming modulations to below 8 Hz.
%
%   Any of the optinal parameters for AUDITORYFILTERBANK, IHCENVELOPE
%   and ADAPTLOOP may be specified for this function. They will be passed
%   to the corresponding functions.
%   
%   No overshoot limit is considered.
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, dau1996preproc, dau1997preproc
%
%   References:
%     A. Kohlrausch, D. Pueschel, and H. Alphei. Temporal resolution and
%       modulation analysis in models of the auditory system. 1992.
%
% Author        : Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/03/2015
% Adapted from  : dau1996preproc.m (code by Torsten Dau, Morten L. Jepsen, Peter L. Sondergaard)
% Last update on: 24/03/2015 % Update this date manually
% Last use on   : 24/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ Warning, Strube model not used ----------
time_pause = 0.1;
disp([mfilename '.m: the code of this function is incorrect, because it uses a GM filterbank instead of the Strube one, be aware of this...'])

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
pause(time_pause)

% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop'};

definput.importdefaults={'ihc_dau','adt_dau'};%'adt_puschel'};%'adt_dau'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
outs.out01_filterbank  = outsig; 

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
outs.out02_ihc         = outsig;

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
outs.out03_adaptloop   = outsig;

% Calculate filter coefficients for the 200 ms (approx.eq to 80 Hz)
% modulation lowpass filter.
% This filter places a pole /very/ close to the unit circle.
Tau_LP = 200e-3; % Dau1996 uses 20e-3 (8 Hz)

mlp_a = exp(-(1/Tau_LP)/fs);
mlp_b = 1 - mlp_a;
mlp_a = [1, -mlp_a];

% Apply the low-pass modulation filter.
outsig = filter(mlp_b,mlp_a,outsig);

% Apply final resampling to avoid excessive data
if ~isempty(keyvals.subfs)
  outsig = fftresample(outsig,round(length(outsig)/fs*subfs));
end;
outs.out04_LPF = outsig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
