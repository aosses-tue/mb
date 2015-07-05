function a3xslt(in,xsl,out)
% function a3xslt(in,xsl,out)
%
% xslt6p1(xml_in,xsl_in,html_out)
%
% All three arguments are file names.
%
% The versions of xerces.jar and saxon.jar which ship with MATLAB 6.1
% must be on your classpath
%
% See also: a3localsettings.m
% 
% Created by the APEX3 Team
% Edited by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created in    : 2013-2014
% Last update on: 01/07/2015 % Update this date manually
% Last use on   : 01/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isunix
    if in(1) == '~'
        error('Please specify the input filename considering the full path, i.e., use ''/home/alejandro/Documents/..'' instead of ''~/Documents/..''');
    end
end

% get absolute path
if (in(1)~=delim && in(2)~=':')
    in=which(in);
end
if (xsl(1)~=delim && xsl(2)~=':')
    xsl=which(xsl);
end

if (~exist(in))
    error(['Input xml file not found: ' in]);
end
if (~exist(xsl))
    error('XSLT script not found, check a3localsettings');
end

in = i_url(in);
xsl = i_url(xsl);

tf = javax.xml.transform.TransformerFactory.newInstance;
streamstyle = javax.xml.transform.stream.StreamSource(xsl);
transformer = tf.newTransformer(streamstyle);

streamin = javax.xml.transform.stream.StreamSource(in);
streamout = javax.xml.transform.stream.StreamResult(out);

transformer.setParameter('target', 'parser');

transformer.transform(streamin,streamout);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function url = i_url(filename)

%filename = which(filename);
% This bit might be Windows-specific.
if (isunix)
    extra='';
else
    extra='/';
end
url = [ 'file://' extra strrep(filename,'\','/') ];
