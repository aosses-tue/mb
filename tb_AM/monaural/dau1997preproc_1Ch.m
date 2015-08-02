function [outsig, fc, mfc, extra] = dau1997preproc_1Ch(insig, fs, fc, varargin)
% function [outsig, fc, mfc, extra] = dau1997preproc_1Ch(insig, fs, fc, varargin)
%
% 1. Description:
%       Auditory model from Dau et. al. 1997
%       Usage:  [outsig, fc] = dau1997preproc_1Ch(insig,fs,fc);
%               [outsig, fc] = dau1997preproc_1Ch(insig,fs,fc...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
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
%   1) a gammatone filter bank with 1-ERB spaced filtes.
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
% Author        : Torsten Dau, Morten L. Jepsen, Peter L. Søndergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 07/10/2014 
% Last use on   : 07/10/2014 
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
keyvals.flow    = fc;
keyvals.fhigh   = fc;
[outsig, fc]            = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

extra.insig             = insig;
extra.out01_filterbank  = outsig; 

% 'haircell' envelope extraction
outsig                  = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
extra.out02_ihc         = outsig;

% non-linear adaptation loops (model units, MU)
outsig                  = adaptloop(outsig,fs,'argimport',flags,keyvals);
extra.out03_adaptloop   = outsig;

% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);
extra.out04_modfilterbank = outsig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end