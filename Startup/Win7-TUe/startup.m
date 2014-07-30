% startup.m
%
% Programmed by Alejandro Osses, TU/e 2014
% PC: 		1. TU/e, Win7, MATLAB 2013a
%           2. TU/e, Win7, MATLAB 2009a
% 
% Created on     : 14/05/2014
% Last updated on: 20/05/2014
% Last used on   : 30/06/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bOnline = 0;

if bOnline
    addpath('G:\MATLAB\Utility');
    disp('Using on-line set-up')
else
    addpath('D:\MATLAB-off-line\Utility');
    disp('Using off-line set-up')
end

disp([mfilename '.m: startup file for TU/e MATLAB'])
slCharacterEncoding('UTF-8');
Start_TUe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
