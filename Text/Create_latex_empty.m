function Create_latex_empty(filename)
% function Create_latex_empty(filename)
%
% 1. Description:
%       Creates a simple empty LaTeX file using 'doc_plain.tex' template 
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       Create_latex_empty;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/05/2014
% Last update on: 29/07/2014 % Update this date manually
% Last use on   : 29/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc.lx_Templates = Get_TUe_paths('lx_Templates');
misc.lx_Text      = Get_TUe_paths('lx_Text'     );
templatefile      = 'doc_plain.tex';

type_doc = input('Which kind of document do you want to create? (1 = weekly update): ');

if nargin == 0
    
    switch type_doc
        case 1 % weekly update
            filename = 'update';
        otherwise
            filename = input('Type a name for your new LaTeX document (between brackets, no spaces): ');
    end
        
    lx_title = input('Type the title of the document (otherwise type enter): ');
end

p = Get_date;
% p.dd = '12';
% p.mm = '05';
% p.yyyy = '2014'; 
p.template  = templatefile;
p.style     = templatefile; % style file
p.filename  = filename;
p.title     = lx_title;

switch type_doc
    case 1
        p.comments = 'MATLAB script used: % complete this manually';
    otherwise
        p.comments = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_latex_info = Get_my_latex_info;

p.author    = my_latex_info.author; % p.author    = 'Alejandro Osses';
p.email     = my_latex_info.email;
p.supervisor= my_latex_info.supervisor; % p.supervisor = 'Armin';
%p.titlehead = my_latex_info.titlehead;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputfile = [misc.lx_Templates templatefile];
outputdir = [misc.lx_Text 'lx' p.yyyy '-' p.mm '-' p.dd '-' filename delim];
outputfile = [outputdir filename '.tex'];

Mkdir(outputdir);

p.functionname=filename;

output = readfile_replace(inputfile,p);

fid=fopen(outputfile, 'w');

fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile ' successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])