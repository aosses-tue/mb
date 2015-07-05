function result=a3transform(filename, script)
% function result=a3transform(filename, script)
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Created by: APEX 3 Team
% Edited by: Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Last edited on: 01/07/2015 (only formatting)
% Last use on   : 01/07/2015
% Last update on: 01/07/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=a3localsettings();

% Use external transformer:
if (0)
    cmdline=[r.xalancmd ' -text -nh -param target \''parser\'' -xsl file:' r.apex_xslt_scripts script ' -in "file:' filename '"']; %'" -out "tempfile2388.txt"'];
    [s, result]=system(cmdline);
end

% Use java transformer:

tempfile=tempname;
a3xslt(filename,[r.apex_xslt_scripts script],tempfile);
result=readfile(tempfile);
delete(tempfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
