function r=a3localsettings
% function r=a3localsettings
%
% 1. Description:
%       Returns APEX 3 toolbox settings.
%       Note: use full paths, i.e., use r.apex_path='/home/alejandro/Documenten/apex/';
%                                   instead of r.apex_path='~/Documenten/apex/';
%             the latter will give problems when using other tb_APEX functions
% 2. Stand-alone example:
%       r = a3localsettings;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Created by: APEX 3 Team
% Edited by: Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Last edited on: 06/03/2015
% Last use on   : 06/03/2015
% Last update on: 01/07/2015 (Linux paths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to xalan, only necessary if APEX does not do the XSLT transformation
r.xalancmd='/usr/bin/xalan';

% Main APEX directory:
if ~isunix
    r.apex_path='C:\Program Files\apex\';
else
    r.apex_path='/home/alejandro/Documenten/apex/'; % updated on 01/07/2015
end

% APEX experiment schema
r.apex_schema=[r.apex_path 'schemas/experiment.xsd'];
% APEX XSLT scripts, only necessary if APEX does not do the XSLT transformation
if ~isunix
    r.apex_xslt_scripts=[r.apex_path 'xslt' delim];
else
    r.apex_xslt_scripts=[r.apex_path 'data' delim 'xslt' delim];
end
% Tool to check XML files, not required
r.xml_check_tool='/usr/bin/xmllint';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
