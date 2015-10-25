function [outsig, fc, outs] = dau1996preproc(insig, fs, varargin)
% function [outsig, fc, outs] = dau1996preproc(insig, fs, varargin)
%
% 1. Description:
%       Auditory model from Dau et. al. 1996.
% 
%   Usage: [outsig, fc] = dau1996preproc(insig,fs);
%          [outsig, fc] = dau1996preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   *Warning:* This code is incorrect, the Dau 1996 models uses a transmission-line 
%   model from Strube 1985, while THIS code erroneously uses the gammatone filters. 
%   If/when the Strube model is included in AMToolbox, this function will be fixed.
%
%   DAU1996PREPROC(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz as described in Dau,
%   Puschel and Kohlrausch (1996a).
%  
%   [outsig,fc]=DAU1996PREPROC(...) additionally returns the center 
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
%   to the corresponding functions.
%
%   The model implemented in this file is not identical to the model
%   published in Dau et. al. (1996a). An overshoot limit has been added to
%   the adaptation stage to fix a problem where abrupt changes in the
%   input signal would cause unnaturally big responses. This is described
%   in Dau et. al. (1997a).
%
% 2. Additional information:
%       See also: auditoryfilterbank, ihcenvelope, adaptloop, dau1997preproc
%       AMT Tollbox version: 0.9.5-0.9.7
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
% Examples:
%       [insig,fs] = greasy;
%       insig = resample(insig,44100,fs);
%       fs = 44100;
%       [outsig, fc] = dau1996preproc(insig, fs);
% 
% Author        : Torsten Dau, Morten L. Jepsen, Peter L. Soendergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 15/10/2014 
% Last use on   : 23/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ------ Warning, Strube model not used ----------
% error(['This code of this function is incorrect. Please see the description ' ...
%        'in the help text.']
   
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

definput.importdefaults={'ihc_dau','adt_dau'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

%% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
if nargout == 3
    outs.out_filterbank  = outsig; 
end

%%% Added by AO, up to L124
do_internal_noise = 0;
if do_internal_noise
    for i = 1:size(outsig,2)
        if i==1
            warning('temporal arrangement of internal noise')
        end
        
        dur = size(outsig,1)/fs;
        n = AM_random_noise(0,fs/2,9.7+3,dur,fs); % 9.7 dB SPL per band
        outsig(:,i) = outsig(:,i)+n;
    end
end
%%%

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
if nargout == 3
    extra.out_ihc = outsig;
end

%% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
if nargout == 3
    extra.out_adaptloop   = outsig;
end

%% Modulation low-pass filter
% Calculate filter coefficients for the 20 ms (approx.eq to 8 Hz) modulation
% lowpass filter.
% This filter places a pole /very/ close to the unit circle.
mlp_a = exp(-(1/0.02)/fs);
mlp_b = 1 - mlp_a;
mlp_a = [1, -mlp_a];

% Apply the low-pass modulation filter.
outsig = filter(mlp_b,mlp_a,outsig);

% Apply final resampling to avoid excessive data
if ~isempty(keyvals.subfs)
  outsig = fftresample(outsig,round(length(outsig)/fs*subfs));
end

if nargout == 3
    extra.out_LPF = outsig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
