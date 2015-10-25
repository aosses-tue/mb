function [outsig, fc, mfc, outs] = jepsen2008preproc(insig, fs, varargin)
% function [outsig, fc, mfc, outs] = jepsen2008preproc(insig, fs, varargin)
%
% 1. Description:
%   Auditory model from Jepsen et. al. 2008
%   Usage: [outsig, fc, mfc, IntRep] = jepsen2008preproc(insig,fs);
%          [outsig, fc, mfc, IntRep] = jepsen2008preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
% 
%   JEPSEN2008PREPROC(insig,fs) computes the internal representation of the 
%   signal insig sampled with a frequency of fs Hz as described in Jepsen, 
%   Ewert and Dau (2008).
%  
%   [outsig,fc]=jepsen2008(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The Jepsen 2008 model consists of the following stages:
% 
%     1) A heaphone filter to simulate the effect of a standard set of
%        headphones (drnl_CASP_debug.m).
%
%     2) A middle ear filter to simulate the effect of the middle ear, and
%        to convert to stapes movement (drnl_CASP_debug.m).
%
%     3) DRNL - Dual resonance non-linear filterbank (drnl_CASP_debug.m).
%
%     4) An envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz (this code, stage 3).
%
%     5) An expansion stage (this code, stage 4).
%
%     6) An adaptation stage modelling nerve adaptation by a cascade of 5
%        loops (adaptloop.m).
%
%     7) A modulation filterbank (modfilterbank.m) or low-pass modulation
%        filter. For the latter type, use 'lowpass' as input argument for 
%        this function.
%
%   Any of the optinal parameters for DRNL, IHCENVELOPE and
%   ADAPTLOOP may be optionally specified for this function. They will be
%   passed to the corresponding functions.
%
% 2. Additional info:
%   See also: drnl_CASP, drnl, ihcenvelope, adaptloop, modfilterbank, dau1997preproc
%
%   References:
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am.,
%     124(1):422-438, 2008.
%     
%   Url: http://amtoolbox.sourceforge.net/doc/monaural/jepsen2008preproc.php
% 
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak. 
% This file is part of AMToolbox version 0.9.5
%
% 3. Stand-alone example:
%       % Examples:
%       [insig,fs] = greasy;
%       insig = resample(insig,44100,fs);
%       fs = 44100;
%       [outsig, fc, mfc] = jepsen2008preproc(insig, fs);
% 
% Author        : Torsten Dau, Morten Loeve Jepsen, Peter L. Soendergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 16/08/2015 
% Last use on   : 21/08/2015 
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

definput.import        ={'drnl_CASP' ,'ihcenvelope','adaptloop'};
definput.importdefaults={'jepsen2008','ihc_jepsen' ,'adt_jepsen'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

%% 1. Up- or down-sampling to fs = 44100 Hz:
% Note copied from Two!Ears code: up- or down-sampling because outer-middle
%       ear filter coefficients are calculated (optimised) for input signals
%       at 44100 Hz
if fs < 44100    
    insig = resample(insig,44100,fs);
    fs = 44100;
elseif fs > 44100
    resampleLen = floor(length(insig)) / (44100/fs);
    insig = resample(insig,fs,44100);
    insig = insig(1:resampleLen);
end

%% 2. DRNL and compensation for middle-ear (middle ear)
[outsig, fc, params] = drnl_CASP_debug(insig, fs, 'argimport',flags,keyvals);

if nargout >= 4
    outs.out_filterbank = outsig;
end

%% 3. 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

%% 4. Gain + Expansion stage
outsig = gaindb(outsig,keyvals.gain_after_drnl); % default AMT is 50 dB
outsig = outsig.^2;

%% 5. non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

if flags.do_absolutethreshold
    warning('absolute threshold only implemented in 1Ch version')
end

%% 6. Downsampling (of the internal representations)
if flags.do_resample_intrep
    % In case of downsampling:
    resampleFac = 4;
    resampleLen = floor(length(insig) / resampleFac);
    fs_intrep = fs / resampleFac;
    
    outsig = resample(outsig,1,resampleFac);
    outsig = outsig(1:resampleLen,:);
else
    % In case of no-resampling:
    fs_intrep = fs;
end

%% 7. Modulation filterbank or modulation low-pass filter
if flags.do_modfilterbank
    % Modulation filterbank
    [outsig,mfc] = modfilterbank(outsig,fs_intrep,fc);
end

if flags.do_lowpass
    % Low-pass modulation filter.
    timeconstant = 20e-3;
    f0 = 1/(2*pi*timeconstant);
    [mlp_b mlp_a] = IRIfolp(f0,fs_intrep);
    mfc = f0;
    outsig = filter(mlp_b,mlp_a,outsig);
end

if nargout >= 4
    outs.fs_intrep = fs_intrep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
