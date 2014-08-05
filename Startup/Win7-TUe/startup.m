% startup.m
%
% Programmed by Alejandro Osses, TU/e 2014
% PC: 		1. TU/e, Win7, MATLAB 2013a
%           2. TU/e, Win7, MATLAB 2009a
% 
% Create on     : 14/05/2014
% Last update on: 28/07/2014 
% Last use on   : 28/07/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bOnline = 0;

if bOnline
    addpath('G:\MATLAB\Utility');
    disp('Using on-line set-up')
else
    % addpath('D:\MATLAB-off-line\Utility');
    addpath('D:\MATLAB_git\Utility');
    disp('Using off-line set-up')
end

disp([mfilename '.m: startup file for TU/e MATLAB'])
slCharacterEncoding('UTF-8');
Start_TUe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%