function result=a3header
% function result=a3header
%
% 1. Description:
%       Creates the header for an APEX experiment (v.3.0).
% 
% 2. Stand-alone example:
%       result = a3header;
%       disp(result)
% 
% Comments by Alejandro Osses V.
% Last edited on: 25/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isunix)
t.encoding='ISO-8859-1';
else
t.encoding='Windows-1252';
end

template = 'header.xml';
result=readfile_replace(template,t);

if length(result) == 0
    template = ['a3templates' delim 'header.xml'];
    result=readfile_replace(template,t);
end