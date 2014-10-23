function Create_latex_figures(p)
% function Create_latex_figures(p)
%
% 1. Description:
%       Creates a LaTeX file adding all the *.eps figures of the selected
%       folder
%       Conditions tested successfully (amount of figures): 1
% 
% 2. Additional info:
%   2.1 Tested cross-platform: No
%   2.2 Dependencies:
%       - Get_TUe_paths, specifying the following valid directories:
%           - 'lx_Templates'
%           - 'lx_Text'
%       - Get_date
%       - Mkdir
%       - latex-fig-template-header.tex file (in lx_Text)
%       - latex-fig-template-footer.tex
%       - tb_APEX (whole directory)
%       - Ensure_field
% 
% 3. Stand-alone example:
%       p.FigWidth = 1;
%       Create_latex_figures(p); % when asked, select the figures folder
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 16/05/2014
% Last update on: 30/06/2014 % Update this date manually
% Last used on  : 22/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    p = [];
end

disp('Select directory where *.eps files are located...')
figurepath      = uigetdir;
figurepath      = [figurepath delim];
suffix          = 'eps';

misc.lx_Templates = Get_TUe_paths('lx_Templates');
misc.lx_Text      = Get_TUe_paths('lx_Text'     );
% style             = 'wu.tex';
style             = 'doc_plain.tex';
templateheader    = 'latex-fig-template-header.tex';
templatefooter    = 'latex-fig-template-footer.tex';

headerfile  =   readfile([misc.lx_Templates 'latex-fig-template-header.tex']);
footerfile  =   readfile([misc.lx_Templates 'latex-fig-template-footer.tex']);

filename = input('Type a name for your new LaTeX document (between brackets, no spaces): ');
lx_title = input('Type the title of the document (otherwise type enter): ');

pt = Get_date;

p.dd    = pt.dd;
p.mm    = pt.mm;
p.yyyy  = pt.yyyy; 
p.templateheader = templateheader;
p.templatefooter = templatefooter;
p.style     = style;
p.filename  = filename;
p.subject   = '';
p.title     = lx_title;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_latex_info = Get_my_latex_info;

p.author    = my_latex_info.author;
p.email     = my_latex_info.email;
p.supervisor = my_latex_info.supervisor;
% p.titlehead = my_latex_info.titlehead;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputfile = [misc.lx_Templates templatefile];
outputdir = [misc.lx_Text 'lx' p.yyyy '-' p.mm '-' p.dd '-' filename delim];
outputfile = [outputdir filename '.tex'];

Mkdir( outputdir);
Mkdir([outputdir 'Figures' delim])

p.functionname=filename;

p = Ensure_field(p,'FigWidth'       ,1);

all_fig_files = dir(fullfile(figurepath,'*.eps'));

for i=1:length(all_fig_files)
    copyfile([figurepath all_fig_files(i).name],[outputdir 'Figures' delim all_fig_files(i).name]);
end

latexfile   =   fopen(outputfile,'w');
disp(['File created: ' outputfile])

fwrite(latexfile,headerfile);
fprintf(latexfile,'\n\\section{Figures}\n\n');
fprintf(latexfile,'Listed files:\\vspace{12pt}\n\n');

for i = 1:length(all_fig_files)
    fprintf(latexfile,['List of figures: %s\n\n'], [name2figname(all_fig_files(i).name(1:end-length(suffix))),'\vspace{3pt}']);
end
fprintf(latexfile,'\n\n \\newpage');

for i = 1:length(all_fig_files)
    
    fprintf(latexfile,'\n\\begin{figure}\n');
    fprintf(latexfile,'\t\\centering\n');
        
    fprintf(latexfile,['\t\\includegraphics[width=' num2str(p.FigWidth) '\\textwidth]{%s}\n'], ['Figures/'  all_fig_files(i).name]);
    
    fprintf(latexfile,'\\end{figure}\n');
    % fprintf(latexfile,'\\newpage\n\n');
    
end

fwrite(latexfile,footerfile);
fclose(latexfile);

% inputfile = latexfile;
output = readfile_replace(outputfile,p);

fid=fopen(outputfile, 'w');

fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile '.m successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end