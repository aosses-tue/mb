function result=a3header

if (isunix)
t.encoding='ISO-8859-1';
else
t.encoding='Windows-1252';
end


result=readfile_replace('header.xml',t);