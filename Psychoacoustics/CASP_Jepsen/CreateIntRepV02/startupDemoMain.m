%--------------------------------------------------------------------------
% startup.m
%--------------------------------------------------------------------------
%   setting the internal MATLAB path 
%
% Author        : C.T. Iben
% Created on    : 20/01/2013
% Received on   : 14/10/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 15/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    Current_dir = fileparts(which(mfilename));
    addpath(genpath(Current_dir));
catch
    addpath(genpath(cd));
end

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eof