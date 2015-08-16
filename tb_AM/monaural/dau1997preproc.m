function [outsig, fc, mfc] = dau1997preproc(insig, fs, varargin)
% function [outsig, fc, mfc] = dau1997preproc(insig, fs, varargin)
%
% 1. Description:
%       Auditory model from Dau et. al. 1997
%       Usage:  [outsig, fc] = dau1997preproc(insig,fs);
%               [outsig, fc] = dau1997preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal (column vector).
%     fs     : sampling rate [Hz].
%  
%   Output parameters:
%     output : output of the modulation filterbank. It is a cell array with 
%              dimensions 31 x 1, so output{i} corresponds to the output 
%              of the modulation filter centred at fc(i). 
% 
%   DAU1997PREPROC(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz as described in Dau,
%   Puschel and Kohlrausch (1997a).
%  
%   [outsig,fc,mfc]=DAU1997PREPROC(...) additionally returns the center
%   frequencies of the filter bank and the center frequencies of the
%   modulation filterbank.
%  
%   The Dau 1997 model consists of the following stages:
%   
%   1) a gammatone filter bank with 1-erb spaced filtes.
%
%   2) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 1000 Hz.
%
%   3) an adaptation stage modelling nerve adaptation by a cascade of 5 loops.
%
%   4) a modulation filterbank
%
%   Any of the optinal parameters for AUDITORYFILTERBANK,
%   IHCENVELOPE and ADAPTLOOP may be optionally specified for this
%   function. They will be passed to the corresponding functions.
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, modfilterbank
%
%   Url: http://amtoolbox.sourceforge.net/doc/monaural/dau1997preproc.php
%
%   References: dau1997mapI dau1997mapII
%
% Author        : Torsten Dau, Morten L. Jepsen, Peter L. Soendergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 25/07/2015 % Update this date manually
% Last use on   : 25/07/2015 % Update this date manually
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
definput.importdefaults={'ihc_dau','adt_dau'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank: 
%   insig  = N x  1
%   outsig = N x 31 (31 gammatone filters)
[outsig, fc]    = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% 'haircell' envelope extraction
outsig          = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

% non-linear adaptation loops (model units, MU)
outsig          = adaptloop(outsig,fs,'argimport',flags,keyvals);

% Modulation filterbank
[outsig,mfc]    = modfilterbank1997(outsig,fs,fc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end