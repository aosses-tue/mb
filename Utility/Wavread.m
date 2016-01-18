function [y,Fs,bits,opt_ck] = Wavread(filename)
% function [y,Fs,bits,opt_ck] = Wavread(filename)
%
% 1. Description:
%       Loads an audio file using audioread (available in MATLAB versions
%       released after 2013). Otherwise use wavread.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 23/06/2014
% Last update on: 07/01/2016 
% Last use on   : 07/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mv = Matlab_version;

try
    if mv >= 8.4 % Since MATLAB 2014b
        [y,Fs] = audioread(filename);
        opt_ck.fmt = audioinfo(filename);
        bits = opt_ck.fmt.BitsPerSample;
    else   
        [y,Fs,bits,opt_ck] = audioread(filename);
    end
catch
    [y,Fs,bits,opt_ck] = wavread(filename);
end
disp([mfilename '.m: ' filename ' read'])

if nargout == 0
    try
        fprintf('%s.m: RMS of %.2f [dBFS]\n',mfilename,rmsdb(y))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
