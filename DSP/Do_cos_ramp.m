function outsig = Do_cos_ramp(insig,fs,attack_ms,release_ms)
% function outsig = Do_cos_ramp(insig,fs,attack,release)
%
% 1. Description:
%       Applies a cosine ramp with attack and release times given in [ms]
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 20/08/2014
% Last update on: 20/08/2014 
% Last use on   : 24/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    fs = 44100;
    warning('Using default sampling frequency of 44100 Hz');
end

if nargin < 3
    attack_ms = 25;
end

if nargin < 4
    release_ms = attack_ms;
end

sig_len = length(insig);
r =  cos_ramp(sig_len,fs,attack_ms,release_ms);
try
    outsig = transpose(r).*insig;
catch
    outsig = r.*insig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
