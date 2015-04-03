function [outsig, fc] = auditoryfilterbank(insig, fs, varargin)
% function [outsig, fc] = auditoryfilterbank(insig, fs, varargin)
% 
%   1. Description:
%     	Linear auditory filterbank
% 
%       Usage: [outsig, fc] = auditoryfilterbank(insig,fs);
%              [outsig, fc] = auditoryfilterbank(insig,fs,...);
%
%       Input parameters:
%           insig  : input acoustic signal.
%           fs     : sampling rate.
% 
%       Output parameters:
%           fc     : centre frequencies
%       outsig : output filtered signal
%  
%   AUDITORYFILTERBANK(insig,fs) applies an auditory filterbank to the
%   input signal insig sampled with a frequency of fs Hz. The filterbank
%   is composed of gammatone filters with 1 ERB wide filters.
%  
%   [outsig,fc]=AUDITORYFILTERBANK(...) additionally returns the center
%   frequencies of the filter bank.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%     'basef',basef  Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value of [] means
%                    no default.
%
%     'langendijk'   Use rectangular filters as in Langendijk (2002).        
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/auditoryfilterbank.php
% 
% Author        : Peter L. S�ndergaard
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 18/03/2014 % Update this date manually
% Last use on   : 02/04/2015 % Update this date manually
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

definput.import={'auditoryfilterbank'}; % Default input

[flags,keyvals,flow,fhigh]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, keyvals.bwmul, keyvals.basef);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs, 'complex');

% Apply the Gammatone filterbank
outsig = 2*real(ufilterbankz(gt_b,gt_a,insig));