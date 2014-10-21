function Create_latex_empty(filename)
% function Create_latex_empty(filename)
%
% 1. Description:
%       Creates a simple empty LaTeX file using 'doc_plain.tex' template. 
%       This template uses the format defined in:
%           - 'style_doc_plain.tex'     if it is not a weekly update
%           - 'style_doc_plain_wu.tex'  if it is a weekly update
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       Create_latex_empty;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/05/2014
% Last update on: 17/07/2014 % Update this date manually
% Last use on   : 17/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc.lx_Templates = Get_TUe_paths('lx_Templates');
misc.lx_Text      = Get_TUe_paths('lx_Text'     );

type_doc = input('Which kind of document do you want to create? (1 = weekly update): ');

if nargin == 0
    
    switch type_doc
        case 1 % weekly update
            filename = 'update';
            templatefile = 'doc_plain.tex';
            style        = 'doc_plain_wu.tex';
        otherwise
            filename = input('Type a name for your new LaTeX document (between brackets, no spaces): ');
            templatefile = 'doc_plain.tex';
            style        = templatefile;
    end
        
    lx_title = input('Type the title of the document (otherwise type enter): ');
end

p = Get_date;
% p.dd = '12';
% p.mm = '05';
% p.yyyy = '2014'; 
switch type_doc
    case 1
        p.subject = 'Weekly update';
    otherwise
        p.subject = '';
end

p.template  = templatefile;
p.style     = style; % style file
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
p.university = my_latex_info.university;
p.department = my_latex_info.department;
p.researchgroup = my_latex_info.researchgroup; 
p.address = my_latex_info.address;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputfile = [misc.lx_Templates templatefile];
outputdir = [misc.lx_Text 'lx' p.yyyy '-' p.mm '-' p.dd '-' filename delim];
outputfile = [outputdir filename '.tex'];

status = Mkdir(outputdir);

if status == 0
    error('Directory already exists')
end

p.functionname=filename;

output = readfile_replace(inputfile,p);

fid=fopen(outputfile, 'w');

fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile ' successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])