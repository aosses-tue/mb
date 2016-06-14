function [score trialID pres_order] = parseWAE_ABC(resultfile)
% function [score trialID pres_order] = parseWAE_ABC(resultfile)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    resultfile = 'D:\Documenten-TUe\06-Lectures-TUe\BEP-2015-2016\03-Odin\Experiment\webaudioevaluationtool-45df85a336d4\saves\save-AO-proc.xml';
    warning('Using default input file...')
end

theStruct = parseXML(resultfile);

count = 0;
% Extract trials (or 'page')
for i = 1:length(theStruct.Children)
    if strcmp(theStruct.Children(i).Name,'page')
        count = count + 1;
        idx(count) = i;
    end
end

count = 0;
for i = idx
    count = count+1;
    thisPage = theStruct.Children(i);
    trialID{count} = thisPage.Attributes(2).Value; % Attr1 = 0; Attr2 = 1; Attr3 = 2;
    
    countthis = 1;
    for k = 1:length(thisPage.Children)
        % thisPage.Children(1).Name % metric (NOT USED)
        if strcmp(thisPage.Children(k).Name,'audioelement');
            expr = sprintf('audio%.0f = thisPage.Children(k);',countthis);
            eval(expr);
            Childnr(countthis) = k;
            countthis = countthis + 1; 
        end
    end
    
    %%%
    for k = 1:length(audio1.Children)
        if strcmp(audio1.Children(1,k).Name,'value')
            scoretmp(1) = str2num(audio1.Children(1,k).Children.Data);
        end
    end
    for k = 1:length(audio1.Attributes)
        if strcmp(audio1.Attributes(1,k).Name,'presentedId')
            stimorder{count,1} = audio1.Attributes(1,k).Value;
        end
    end
    
    %%%
    for k = 1:length(audio2.Children)
        if strcmp(audio2.Children(1,k).Name,'value')
            scoretmp(2) = str2num(audio2.Children(1,k).Children.Data);
        end
    end
    for k = 1:length(audio2.Attributes)
        if strcmp(audio2.Attributes(1,k).Name,'presentedId')
            stimorder{count,2} = audio2.Attributes(1,k).Value;
        end
    end
    
    %%%
    for k = 1:length(audio3.Children)
        if strcmp(audio3.Children(1,k).Name,'value')
            scoretmp(3) = str2num(audio3.Children(1,k).Children.Data);
        end
    end
    for k = 1:length(audio3.Attributes)
        if strcmp(audio3.Attributes(1,k).Name,'presentedId')
            stimorder{count,3} = audio3.Attributes(1,k).Value;
        end
    end
    
    [xx idx] = sort(stimorder(count,:)); %,'descend');
    score(count,:) = scoretmp(idx);
    
end

pres_order = 1:length(trialID); 

[xx idx] = sort(trialID);
score = score(idx,:);
pres_order = pres_order(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end