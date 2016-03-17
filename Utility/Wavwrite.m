function Wavwrite(x,fs,outputfilename, varargin)
% function Wavwrite(x,fs,outputfilename, varargin)
% 
% 1. Description:
%       Write Microsoft WAVE (".wav") sound file
% 
% 2. Stand-alone example:
%       % Use it the same way as the old 'wavwrite' function.
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 13/05/2014
% Last update on: 18/06/2015 
% Last use on   : 16/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <= 3
    try
        audiowrite(outputfilename,x,fs);
    catch
        wavwrite(x,fs,outputfilename);
    end
elseif nargin == 4
    nbits = outputfilename;
    outputfilename = varargin{1};
    try
        audiowrite(outputfilename,x,fs,'BitsPerSample',nbits);
    catch
        wavwrite(x,fs,nbits,filename);
    end
end

disp([mfilename '.m: file ' outputfilename ' created'])
try
    fprintf('%s.m: RMS of %.2f [dBFS]\n',mfilename,rmsdb(x))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
