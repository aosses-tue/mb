function [y] = Zero_padding(x,time2pad,fs)
% function [y] = Zero_padding(x,time2pad,fs)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/08/2014
% Last update on: 16/08/2014 % Update this date manually
% Last use on   : 16/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    fs=44100;
end    

if nargin < 2
    time2pad = 200e-3; % 200 ms
end    

N = round(fs*time2pad);

if size(x,1) == 1
    y = [zeros(1,N) x zeros(1,N)]; 
end

if size(x,2) == 1
    y = [zeros(N,1); x]; % y = [zeros(N,1); x; zeros(N,1)]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
