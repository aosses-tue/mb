function [outsig, fc] = dorp2011preproc(insig, fs, varargin);
% function [outsig, fc] = dorp2011preproc(insig, fs, varargin);
%
% 1. Description:
%  Auditory model from van Dorp et. al. 2001
%   Usage: [outsig, fc] = dorp2011preproc(insig,fs);
%          [outsig, fc] = dorp2011preproc(insig,fs,...);
%
%   Input parameters:
%        insig  : input acoustic signal.
%        fs     : sampling rate.
%  
%   DORP2011PREPROC(insig,fs,tau,ild) computes the EI-cell
%   representation of the signal insig sampled with a frequency of fs Hz
%   as described in Breebaart (2001). 
%
%   The input must have dimensions time x left/right channel
%   x signal no.
%
%   The output has dimensions time x frequency x signal no. 
%  
%   [outsig,fc]=DORP2011PREPROC(...) additionally returns the center
%   frequencies of the filter bank.
%  
%   The Breebaart 2001 model consists of the following stages:
%   
%   1) a gammatone filter bank with 1-erb spaced filters.
%
%   2) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz.
%
%   3) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops.
%
%   4) an excitation-inhibition (EI) cell model.
%
%   Parameters for AUDITORYFILTERBANK, IHCENVELOPE, ADAPTLOOP and
%   EICELL can be passed at the end of the line of input arguments.
%
%   Examples
%   --------
%
%   The following code sets up a simple test example :
%
%     % Setup parameters
%     fs      = 44100;            % Sampling rate
%     T       = 0.3;              % Duration
%     Spl1    = 75;               % SPL of input signal 1
%     Spl2    = 75;               % SPL of input signal 2
%     rho     = 0;                % normalized correlation of signals
%     tau     = 0;
%     ild     = 0;
%
%     % Generate signals:
%     t  = [0:1/fs:T];
%     n1 = setdbspl(randn(length(t),1),Spl1);
%     n2 = setdbspl(randn(length(t),1),Spl2);
%     x1 = n1*sqrt((1+rho)/2) + n2*sqrt((1-rho)/2);
%     x2 = n1*sqrt((1+rho)/2) - n2*sqrt((1-rho)/2);
%
%     % % Run the model and plot it
%     [ei_map, fc] = dorp2011preproc([x1,x2], fs);
%     % plotfilterbank(ei_map,1,fc,fs,'audtick','lin');
%       file = 'D:\Databases\dir01-Instruments\Piano\00-Original-files\pressionexpeCd5.wav';
%       [insig fs] = Wavread(file);
%       lvl = rmsdb(insig)+100;
%       [par psi] = raa(file,100);
% 
%   See also: eicell, auditoryfilterbank, ihcenvelope, adaptloop
%
%   Url: http://amtoolbox.sourceforge.net/doc/binaural/breebaart2001preproc.php

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
% 
%   References: breebaart2001binaural

%   AUTHOR : Peter L. Soendergaard
  
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

definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop'};
definput.importdefaults={'gtf_dorp','ihc_breebaart','adt_dorp'};

[flags,keyvals,flow,fhigh,basef]  = ltfatarghelper({'flow','fhigh','basef'},definput,varargin);

% ------ do the computation -------------------------

%% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

%% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

% [siglen,nfreqchannels,naudiochannels,nsignals] = size(outsig);

[mlp_b mlp_a] = IRIfolp(8,fs);
outsig = filter(mlp_b,mlp_a,outsig);
    
if nargout == 0
    
    idx_i = 5-4;
    idx_f = 9-4;
    K = (idx_f-idx_i)+1;
    N = size(outsig,1);
    LLow = 1/(K*N)*sum(  sum( sqrt(outsig(:,idx_i:idx_f,1).^2+outsig(:,idx_i:idx_f,2).^2 ))  ); % Eq. 3.25, pp 82
    
    disp('Energy per channel and per band...')
    mean( abs(outsig) )
end