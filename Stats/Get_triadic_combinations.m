function [triads Ncombinations pairs] = Get_triadic_combinations(minM,maxM)
% function [triads Ncombinations pairs] = Get_triadic_combinations(minM,maxM)
%
% 1. Description:
%       Get pairwise combinations for an array going from minM to maxM in
%       steps of 1.
% 
% 2. Stand-alone example:
%       minM = 1;
%       maxM = 5;
%       Get_triadic_combinations(minM,maxM);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 27/04/2016
% Last update on: 27/04/2016 
% Last use on   : 27/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ms      = minM:maxM;
N       = length(Ms);
triads   = [];
for i = 1:N
    for j = 1:N
        for k = 1:N
            triads(end+1,:) = sort([i j k]);
        end
    end
end

for i = size(triads,1):-1:1
    if triads(i,1) == triads(i,2) | triads(i,2) == triads(i,3)
        triads(i,:) = [];
    end
end

triads = unique(triads,'rows');
Ncombinations = size(triads,1); % nchoosek(maxM - minM+1,3)

if nargout == 3
    pairs = Get_pairwise_combinations(minM,maxM);
    countpair = 0;
    for i = 1:size(triads,1)
        if ( sum(pairs(1,1)==triads(i,:)) + sum(pairs(1,2)==triads(i,:)) )== 2
            countpair = countpair + 1;
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
