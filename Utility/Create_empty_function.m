function Create_empty_function(filename, options)
% function Create_empty_function(filename, options)
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/05/2014
% Last update on: 02/08/2014 % Update this date manually
% Last used on  : 02/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    options = [];
end

options = Ensure_field(options,'bIncludeLog',input('Do you want to include a log for the script/function? (1 = yes, 0 = no): '));

misc.MATLAB = Get_TUe_paths('MATLAB');

if nargin == 0
    filename = input('Type a name for your new function (between brackets, no spaces): ');
end

inputfile = [misc.MATLAB 'template_function.m'];
outputfile = [misc.MATLAB 'Utility' delim 'New-scripts' delim filename '.m'];
Mkdir([misc.MATLAB 'Utility' delim 'New-scripts']);

p = Get_date;

p.functionname=filename;

if options.bIncludeLog
    p.log_begin = sprintf('\nbDiary = 1;\nDiary(mfilename,bDiary);\n\n');
    p.log_end   = sprintf('\nif bDiary\n\tdiary off\nend\n');
    p.eof       = 'disp([''EOF: '' mfilename ''.m''])';
else
    p.log_begin = '';
    p.log_end   = '';
    p.eof       = 'end';
end

output = readfile_replace(inputfile,p);

fid=fopen(outputfile, 'w');

fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile '.m successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
