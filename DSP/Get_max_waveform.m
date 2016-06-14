function Nmax = Get_max_waveform(insig,fs,winlen)
% function Nmax = Get_max_waveform(insig,winlen)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 07/05/2016
% Last update on: 07/05/2016 
% Last use on   : 07/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamples = 1:length(insig);
t        = Nsamples/fs;
[RMS_se1 t_se1] = rmsdb_sec(insig,fs,winlen,0);

%   2.1. Detecting the maximum and matching the levels:
[max_1 idx] = max(RMS_se1);
idx = find(t < t_se1(idx),1,'last');

Nmax = Nsamples(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
