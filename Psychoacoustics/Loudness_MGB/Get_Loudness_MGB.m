function [InstantaneousLoudness, ShortTermLoudness, LongTermLoudness, times] = Get_Loudness_MGB(insig, fs, filtermethod)
% function [InstantaneousLoudness, ShortTermLoudness, LongTermLoudness, times] = Get_Loudness_MGB(insig, fs, filtermethod)
%
% 1. Description:
%
% 2. Stand-alone example:
%       filename = 'D:\MATLAB\NS19-Cd5_2.wav'; 
%       [insig fs] = Wavread(filename);
%       Get_Loudness_MGB(insig, fs);       
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 29/03/2016
% Last update on: 29/03/2016 
% Last use on   : 29/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    filtermethod = 1;
end

if nargout == 0
    doplot = 1;
else
    doplot = 0;
end
cal = 34.0507; % set to make 1 kHz tone at 40 dB yield 1 sone
insig = From_dB(6)*insig; % if cal = 34.0507

faster = 1; % if you think this analyser is slow, try setting this to 0!
decay = 5000; % 5 seconds of decay after the end of the file

[InstantaneousLoudness, ShortTermLoudness, LongTermLoudness, times] = MGBLoudness2b(insig,fs,filtermethod,cal,faster,decay,doplot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
