function outsig = Randomise_insig(insig)
% function outsig = Randomise_insig(insig)
%
% 1. Description:
%       outsig containes the same samples (same signal length) than insig1 
%       but with a starting sample that is randomly defined (following a 
%       uniform distribution). 
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/08/2015
% Last update on: 12/08/2015 
% Last use on   : 12/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
