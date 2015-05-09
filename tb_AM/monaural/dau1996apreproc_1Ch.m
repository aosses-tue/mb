function [outsig, fc, outs] = dau1996apreproc_1Ch(insig, fs, fc, varargin);
% function [outsig, fc, outs] = dau1996apreproc_1Ch(insig, fs, fc, varargin);
%
% 1. Description:
%       Auditory model from Dau et. al. 1996
% 
%   Usage: [outsig, fc, outs] = dau1996apreproc(insig,fs,fc);
%          [outsig, fc, outs] = dau1996apreproc(insig,fs,fc,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%     outs   : struct with outputs of each stage of the model (added by AO)
%
%   *Warning:* This code is incorrect, the Dau 1996 models uses a transmission-line 
%   model from Strube 1985, while THIS code erroneously uses the gammatone filters. 
%   If/when the Strube model is included in AMToolbox, this function will be fixed.
%
%   DAU1996APREPROC(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz as described in Dau,
%   Puschel and Kohlrausch (1996a).
%  
%   [outsig,fc]=DAU1996APREPROC(...) additionally returns the center 
%   frequencies of the filter bank.
%
%   The Dau 1996 model consists of the following stages:
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
%   to the corresponding functions. No overshoot limit is considered.
%   
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, dau1996preproc, dau1997preproc
%
%   References:
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615-3622, 1996a.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     "effective" signal processing in the auditory system. II. Simulations
%     and measurements. J. Acoust. Soc. Am., 99:3623-3631, 1996b.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892-2905, 1997a.
%
%   Url: http://amtoolbox.sourceforge.net/doc/monaural/dau1996preproc.php
%
% Author        : Torsten Dau, Morten L. Jepsen, Peter L. Sondergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 08/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

definput.importdefaults={'ihc_dau','adt_dau1996'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
keyvals.flow    = fc;
keyvals.fhigh   = fc;
[outsig, fc]            = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
outs.out01_filterbank   = outsig; 

% 'haircell' envelope extraction
outsig                  = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
outs.out02_ihc          = outsig;

% non-linear adaptation loops
outsig                  = adaptloop(outsig,fs,'argimport',flags,keyvals);
outs.out03_adaptloop    = outsig;

% Calculate filter coefficients for the 20 ms (approx.eq to 8 Hz) modulation
% lowpass filter.
% This filter places a pole /very/ close to the unit circle.
mlp_a   = exp(-(1/0.02)/fs);
mlp_b   = 1 - mlp_a;
mlp_a   = [1, -mlp_a];

% Apply the low-pass modulation filter.
outsig = filter(mlp_b,mlp_a,outsig);

% Apply final resampling to avoid excessive data
if ~isempty(keyvals.subfs)
  outsig = fftresample(outsig,round(length(outsig)/fs*subfs));
end;
outs.out04_LPF = outsig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
