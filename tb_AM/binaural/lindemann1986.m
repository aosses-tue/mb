function [crosscorr,t,ild,cfreq,params] = lindemann1986(insig,fs,varargin)
% function [crosscorr,t,ild,cfreq,params] = lindemann1986(insig,fs,varargin)
%   
%   1. Description:
%   Calculates a binaural activation pattern
% 
%   Usage: [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f,T_int,N_1)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f,T_int)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s)
%          [crosscorr,t] = lindemann1986(insig,fs)
%
%   Input parameters:
%       insig       : binaural signal for which the cross-correlation
%                     should be calculated
%       fs          : sampling rate (Hz)
%       c_s
%       w_f
%       M_f         : fading constant for the monaural sensitivity (default = 6)
%       T_int
%       N_1         : lower summation boundary for the running 
%                     cross-correlation (default = 1)
%
%   Output parameters:
%       crosscorr   : A matrix containing the cross-correlation signal
%                     for every frequency channel fc and every time step n.
%                     The format of this matrix is output(n,m,fc), where m*
%                     denotes the correlation (delay line) time step.
%       t           : time axis for the time steps n in crosscorr
%       ild         : interaural level difference (ILD) for every freqeuncy
%                     channel fc*
%       cfreq       : center frequencies of every frequency channel
%       params      : struct with input parameters. If any input param was 
%                     omitted then the default value will be returned
%
%   LINDEMANN1986(insig,fs) calculates a binaural activity map for the given
%   insig using a cross-correlation (delay-line) mechanism. The calculation
%   is done for every frequency band in the range 5-40 Erb.
%
%   Lindemann has extended the delay line model of Jeffres (1948) by a
%   contralateral inhibition, which introduce the ILD to the model.  Also
%   monaural detectors were extended, to handle monaural signals (and some
%   stimuli with a split off of the lateralization image). Hess has
%   extented the output from the Lindemann model to a binaural activity map
%   dependend on time, by using a running cross-correlation function.
%   This has been done here by starting a new running cross-correlation
%   every time step T_int.  A detailed description of these cross-
%   correlation steps is given in the LINDEMANN1986BINCORR function.
%
%   The steps of the binaural model to calculate the result are the
%   following:
%
%   1) The given stimulus is filtered using an erb bank to
%      get 36 frequency bands containing a stimulus waveform.
%
%   2) In a second step the auditory nerve is siumulated by extracting the
%      envelope using a first order low pass filter with a cutoff frequency
%      of 800 Hz and half-wave rectification.
%
%   3) Calculation of the cross-correlation between the left and right
%      channel.  This is done using the model described in Lindemann
%      (1986a) and Hess (2007). These are extensions to the delay line model
%      of Jeffres (1948).
%
%   You may supply any flags or key/value pairs of the AUDITORYFILTERBANK,
%   IHCENVELOPE or LINDEMANN1986BINCORR at the end of the line of input
%   arguments.
%
%   Examples:
%   ---------
%
%   This example shows how to the binaural activity map for one frequency
%   channel of the Lindemann binaural model for a sinusoid with a binaural
%   modulation rate of 2 Hz. :
%     
%     fs = 44100; % Sampling rate    
%     f = 500;    % Frequency of the sinusoid
%     mf = 2;     % Binaural modulation frequency
%
%     % Generate 1~s binaural modulated sinusoid
%     sig = bmsin(f,mf,fs);
%
%     % Model paramter (Note: T_int (ms) should be a multiple of 1000/f == 2)
%     % Begin of the storage of the cross-correlation is set to 1, because we have a
%     % non-stationary signal
%
%     % Calculate binaural cross-correlation
%     [cc,t] = lindemann1986(sig,fs,'T_int',6);
%
%     % Plot frequency channel 11, due to round(freqtoerb(500))==11
%     plotlindemann1986(cc,t,'fc',f);
%
%   See also: lindemann1986bincorr, plotlindemann1986, gammatone, ufilterbankz
%
%   Demos: demo_lindemann1986
%
%   References:
%     W. Gaik. Combined evaluation of interaural time and intensity
%     differences: Psychoacoustic results and computer modeling. J. Acoust.
%     Soc. Am., 94:98-110, 1993.
%     
%     W. Hess. Time-Variant Binaural-Activity Characteristics as Indicator of
%     Auditory Spatial Attributes. PhD thesis, Ruhr-Universitaet Bochum,
%     2007.
%     
%     L. Jeffress. A place theory of sound localization. Journal of
%     comparative and physiological psychology, 41(1):35-39, 1948.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608-1622, 1986.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. II. The law of the first wave front. J.
%     Acoust. Soc. Am., 80:1623-1630, 1986.
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/binaural/lindemann1986.php

% Copyright (C) 2009-2014 Peter L. Sondergaard and Piotr Majdak.
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

%   A first implementation of the Lindemann model in MATLAB was done from
%   Wolfgang Hess and inspired this work.

% AUTHOR: Hagen Wierstorf
% Edited by: Alejandro Osses, TU/e 2014
% Last edited on: 15/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------ Checking of input  parameters ---------------------------------
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;
if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end
if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

% Parse the command line and load default parameters
% For default values see lindemann1986a page 1613
% NOTE: I modified the default value for T_int from 10 to 5.
proc = {'auditoryfilterbank','ihcenvelope','lindemann1986bincorr'};
definput.import=proc;

% Highest and lowest frequency to use for the ERB-filterbank (this gives us
% 36 frequency channels, channel 5-40)
definput.importdefaults               = {'flow',erbtofreq(5),'fhigh',erbtofreq(40),'ihc_lindemann'};
[flags,keyvals,c_s,w_f,M_f,T_int,N_1] = ltfatarghelper({'c_s','w_f','M_f','T_int','N_1'},definput,varargin);

params.c_s = c_s;
params.w_f = w_f;
params.M_f = M_f;
params.T_int = T_int;
params.N_1 = N_1;
params.proc = definput.import;
params.text = sprintf('Model params: \n c_inh = %.1f \n T_int = %i [ms] \n w_f = %.3f \n M_f = %i',c_s,T_int,w_f,M_f);
disp(params.text);

%% ------ Computation -----------------------------------------------------

% - ERB Bank --------------------------------------------------------------
% Apply the auditory filterbank
% NOTE: Lindemann uses a bandpass filterbank after Duifhuis (1972) and
% Blauert and Cobben (1978).
[inoutsig,cfreq] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
% dim(inoutsig) = N x 36 x 2 = N sample times, 36 filters, 2 channels (L and R)

% - ILD -------------------------------------------------------------------
% Calculate the interaural level difference (ILD) for every frequency channel
% NOTE: this was not part of the original Lindemann model
try
    ild = Get_ILD(inoutsig(:,:,1),(inoutsig(:,:,2)));
    disp('Get_ILD was used')
catch
    ild = dbspl(inoutsig(:,:,2))-dbspl(inoutsig(:,:,1)); % 1 x 36 ILDs, one per band
end

% - Cross-correlation computation -----------------------------------------
% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
% ------ Haircell simulation -------
% Half-wave rectification and envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'argimport',flags,keyvals);
% ------ Cross-correlation ------
% Calculate the cross-correlation after Lindemann (1986a).
[crosscorr,t] = lindemann1986bincorr(inoutsig,fs,c_s,w_f,M_f,T_int,N_1);

