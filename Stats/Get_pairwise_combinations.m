function [pairs Ncombinations] = Get_pairwise_combinations(minM,maxM)
% function [pairs Ncombinations] = Get_pairwise_combinations(minM,maxM)
%
% 1. Description:
%       Get pairwise combinations for an array going from minM to maxM in
%       steps of 1.
% 
% 2. Stand-alone example:
%       minM = 1;
%       maxM = 5;
%       Get_pairwise_combinations(minM,maxM);
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
pairs   = [];
for i = 1:N
    for j = 1:N
        pairs(end+1,:) = sort([i j]);
    end
end

for i = size(pairs,1):-1:1
    if pairs(i,1) == pairs(i,2)
        pairs(i,:) = [];
    end
end

pairs = unique(pairs,'rows');
Ncombinations = (N-1)*N/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
