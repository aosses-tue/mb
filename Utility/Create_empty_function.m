function Create_empty_function(filename)
% function Create_empty_function(filename)
%
% 1. Description:
%   Creates an empty MATLAB function (m-file)
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example to create a function named hola.m:
%   Create_empty_function('hola');
% 
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 12/05/2014
% Last update: 14/05/2014 % Update this date manually
% Last used: 17/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc.MATLAB = Get_TUe_paths('MATLAB');

if nargin == 0
    filename = input('Type a name for your new function (between brackets, no spaces): ');
end

inputfile = [misc.MATLAB 'template_function.m'];
outputfile = [misc.MATLAB 'Utility' delim 'New-scripts' delim filename '.m'];
Mkdir([misc.MATLAB 'Utility' delim 'New-scripts']);

p = Get_date;
% p.dd = '12';
% p.mm = '05';
% p.yyyy = '2014'; 

p.functionname=filename;

output = readfile_replace(inputfile,p);

fid=fopen(outputfile, 'w');

fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile '.m successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])