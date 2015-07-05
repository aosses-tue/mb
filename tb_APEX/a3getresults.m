function [results,parameters,general] = a3getresults(filename,script,forcetransform)
% function [results,parameters,general] = a3getresults(filename,script,forcetransform)
%
% 1. Description:
%       Get results from an APEX3 result file transform first if necessary 
%       (check implementation on windows)
%
% 2. Stand-alone examples:
% % 2.1. Example at ExpORL, KUL:
%       filename = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/ci-Jan_Leys/20131007-PR/PR_Ref_104_UW_ACE-LB-ACE-JL.apr';
%       forcetransform = 1;
%       [results,parameters,general] = a3getresults(filename,script,forcetransform);
%
% % 2.2. Generic example (it should be cross-platform):
%       [file1,file2]=uigetfile('*.*','Select an APEX result file'); 
%       filename = [file2 file1];
%       script = 'apexresult.xsl';
%       forcetransform = 1;
%       [results,parameters,general] = a3getresults(filename,script,forcetransform);
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

if (nargin<2)
    script='';
end
if (nargin<3)
    forcetransform=0;
end

if (~exist(filename, 'file'))
    % look in all data paths
    pa=getpaths;
    paths=pa.results;
    for i=1:length(paths)
        newpath=fixpath([paths{i} '/' filename ]);
        if (exist( newpath ,'file'))
            filename=newpath;
        end
    end
end

results=readpartfile(filename,'<processed>', '</processed>');

if (forcetransform || length(results)==0)
    if (length(script)==0)
        script=a3getxsltscript(filename);
    end
    
    if (length(script)==0)
        error(['no XSLT script found for file ' filename]);
    end
    
    text=a3transform(filename, script);
    results=strread(text, '%[^\n]\n');
end


if (nargout>1)
    parameters=struct;
    part=readpartfile(filename,'<parameters>', '</parameters>');
    for l=1:length(part)
        line=part{l};
        names=regexp(line, '<parameter name="(?<name>[^"]*)">(?<value>[^<]*)<', 'names');
        if (~isempty(names))
            parameters=setfield( parameters, escapefieldname_l(names(1).name), names(1).value);
        end
    end
end

if (nargout>2)
    general=struct;
    part=readpartfile(filename,'<general>', '</general>');
    for l=1:length(part)
        line=part{l};
        names=regexp(line, '<(?<name>[^"]*)>(?<value>[^<]*)<', 'names');
        if (~isempty(names))
            general.(names(1).name) =  names(1).value;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function name=escapefieldname_l(name)
name = regexprep(name, ' ', '_');