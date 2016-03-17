function result=a3results(type, matlabscript,subject,showresults, saveprocessedresults, xsltscriptparameters)
% result=a3results(filename)
% filename = name of xslt script to be used

if (nargin<3)
    subject=NaN;
    subjectDefined=0;
else
    subjectDefined=1;
end
if (nargin<2)
    matlabscript='';
end

if (nargin<4)
    showresults=0;
end
if (nargin<5)
    saveprocessedresults=0;
end
if (nargin<6)
    xsltscriptparameters='';
end

if (length(type)<4 ||  (length(type)>=4) && ~strcmp(type(end-3:end), '.xsl'))

    if (strcmp(type, 'mao'))
        script='mao.xsl';
    elseif (strcmp(type, 'opc'))
        script='opc.xsl';
    elseif (strcmp(type, 'cst'))
        script='cst.xsl';
    else
        error('Unsupported type');
    end
else
    script=type;
end


lf=sprintf('\n');
result=['<results>' lf];
result=[result '<xsltscript>' script '</xsltscript>' lf];
if (length(xsltscriptparameters))
    result=[result wraptag('xsltscriptparameters', xsltscriptparameters) lf];
end

result=[result '<showresults>'  bool2xml(showresults) '</showresults>' lf];
result=[result '<saveprocessedresults>'  bool2xml(saveprocessedresults) '</saveprocessedresults>' lf];
if (length(matlabscript))
    result=[result wraptag('matlabscript', matlabscript) lf];
end
if (subjectDefined)
    result=[result wraptag('subject', subject) lf];
end

result=[result '</results>' lf ];
