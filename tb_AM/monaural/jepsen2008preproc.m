function [outsig, fc, mfc] = jepsen2008preproc(insig, fs, varargin);
%JEPSEN2008PREPROC   Auditory model from Jepsen et. al. 2008
%   Usage: [outsig, fc] = jepsen2008preproc(insig,fs);
%          [outsig, fc] = jepsen2008preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   *Warning:* This code cannot be verified. It has not been possible to
%   tell from the desciption in the original paper nor from personal
%   communication with the original authors what the correct parameter set
%   used for the model is. This code is kept here as a reminder of the
%   structure of the model, and may reappear in a future work if a
%   verified parameter set can be established. The status of this piece
%   of code is "not even wrong": http://en.wikipedia.org/wiki/Not_even_wrong.
%
%   JEPSEN2008PREPROC(insig,fs) computes the internal representation of the signal insig
%   sampled with a frequency of fs Hz as described in Jepsen, Ewert and
%   Dau (2008).
%  
%   [outsig,fc]=jepsen2008(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The Jepsen 2008 model consists of the following stages:
% 
%     1) a heaphone filter to simulate the effect of a standard set of
%        headphones.
%
%     2) a middle ear filter to simulate the effect of the middle ear, and
%        to convert to stapes movement.
%
%     3) DRNL - Dual resonance non-linear filterbank.
%
%     4) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     5) an expansion stage
%
%     6) an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     7) a modulation filterbank.
%
%   Any of the optinal parameters for DRNL, IHCENVELOPE and
%   ADAPTLOOP may be optionally specified for this function. They will be
%   passed to the corresponding functions.
%
%   See also: drnl, ihcenvelope, adaptloop, modfilterbank, dau1997preproc
%
%   References:
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am.,
%     124(1):422-438, 2008.
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/monaural/jepsen2008preproc.php
% 
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Examples:
%       [insig,fs] = greasy;
%       insig = resample(insig,44100,fs);
%       fs = 44100;
%       [outsig, fc, mfc] = jepsen2008preproc(insig, fs);
% 
% Author        : Torsten Dau, Morten Loeve Jepsen, Peter L. Soendergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 22/04/2015 % Update this date manually
% Last use on   : 22/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ Checking of input parameters ------------

warning(['This code of this function is incorrect. Please see the description ' ...
       'in the help text.'])

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'drnl','ihcenvelope','adaptloop'};
% definput.importdefaults={'jepsen2008'};
definput.importdefaults={'jepsen2008','ihc_jepsen','adt_dau'}; % Added by AO
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avoid truncation errors due to resampling: COPIED FROM CASP ALGO.
IntRep.resampleFac = 4;
IntRep.resampleLen = floor(length(insig) / IntRep.resampleFac);
IntRep.fs = fs / IntRep.resampleFac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Headphone filter (outer ear)
%       Typical human headphone-to-eardrum gain. 
%       insig is assumed to be in [Pa]

hp_fir = headphonefilter(fs);       % Getting the filter coefficients at fs
outsig = filter(hp_fir,1,insig);    % Applying the FIR filter

%% 2. DRNL and compensation for middle-ear (middle ear)

% [outsig, fc] = drnl(outsig, fs, 'argimport',flags,keyvals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[outsig, fc] = drnl_CASP(outsig, fs, 'argimport',flags,keyvals,IntRep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

outsig = gaindb(outsig,50); % moved by AO from line after DRNL filterbank

%% 4. Expansion stage
outsig = outsig.^2;

%% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

%% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
