function Wavwrite(x,fs,outputfilename, varargin)
% function Wavwrite(x,fs,outputfilename, varargin)
% 
% 1. Description:
%       Write Microsoft WAVE (".wav") sound file
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/05/2014
% Last update on: 18/06/2015 % Update this date manually
% Last use on   : 18/06/2015 % Update this date manually
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
