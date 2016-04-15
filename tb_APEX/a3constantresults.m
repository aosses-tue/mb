function [sequence,stimuli, procID]=a3constantresults(filename,xsltscript,forcetransform)
% function [sequence,stimuli, procID]=a3constantresults(filename,xsltscript,forcetransform)
%
% 1. Description:
%       Return adaptive staircase from APEX 3 results file.
%
% 2. Stand-alone example:
%       filename = '~/Documenten/Meas/Meas/Experiments/Results_XML/ci-Jean-Baptiste_Daumerie/20131016-LT/SPIN-LIST-vrouw_30-LISTvrouw_ltass-results-baseline-JB.apr'; 
%       [sequence,stimuli]=a3adaptiveresults(filename);
% 
% 3. Additional information:
%       Tested cross-platform: yes
% 
% Programmed by ExpORL, KU Leuven, Belgium, 2012-2014
% Comments by Alejandro Osses V.
% Last used on: 30/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
    xsltscript='apexresult.xsl';
    forcetransform=0;
end
if (nargin==2)
    forcetransform=1;
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


%% use new apex version
results = a3getresults(filename,xsltscript,forcetransform);
s       = a3parseresults(results);

if (~isfield(s(1),'procedure'))
    
    sequence=zeros(length(s),1);
    stimuli=cell(length(s),1);
    for i=1:length(s)
        sequence(i)=str2num(s(i).adaptiveparameter);
        stimuli{i}=s(i).stimulus;
    end
    
else        % Multiprocedure
    % Get procedure IDs
%     procid=cell(length(s),1);
%     for i=1:length(s)
%         procid{i} = s(i).procedure;
%     end
%     uproc = sort(unique(procid));
    
    sequences   = struct;
    
    for i=1:length(s)
        
        proc = s(i).procedure;
        procID{i} = proc;
        stimuli.(proc){i} = {s(i).stimulus};
        
        if (~isfield(sequences, proc))
            sequences.(proc) = [];
        end
        if strcmp(s(i).corrector,'true')
            sequences.(proc) = [sequences.(proc) 1];
        else
            sequences.(proc) = [sequences.(proc) 0];
        end
    end

    sequence = sequences;
end

procID = unique(procID);