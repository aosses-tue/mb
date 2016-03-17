function result=delim
% function result=delim
% 
% 1. Description:
%       Returns file separator or 'delimiter' associated to the operating
%       system in use.
%       - For Unix-based system the result is '/'
%       - For Windows-based system the result is '\'
% 
% 2. Stand-alone example:
%       string = cd;
%       string = [string delim 'New-folder' delim]; % Returns a string, assuming
%                               % a folder called 'New-folder' inside the current
%                               % MATLAB directory
% 
% Programmed at ExpORL, KU Leuveb 2012-2013
% Comments by Alejandro Osses V., HTI Group, TU/e
% Last used on: 16/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isunix)
    result='/';
else
    result='\';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%