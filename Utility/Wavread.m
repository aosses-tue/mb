function [y,Fs,bits,opt_ck] = Wavread(filename)
% function [y,Fs,bits,opt_ck] = Wavread(filename)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/06/2014
% Last update on: 29/07/2014 % Update this date manually
% Last use on   : 29/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,Fs,bits,opt_ck] = wavread(filename);
disp([mfilename '.m: ' filename ' read'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
