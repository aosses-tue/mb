function p = sone2phon(s)
% function p = sone2phon(s)
%
% 1. Description:
%
% 2. Stand-alone example:
%       p = [20 45];
%       s = phon2sone(p);
%       p_after = sone2phon(s);
% 
% 3. Additional info:
%       See also: phon2sone
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/01/2016
% Last update on: 04/01/2016 
% Last use on   : 04/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx1 = find(s >= 1);
idx2 = find(s < 1);

p(idx1) = 10*log10(s(idx1))/log10(2) + 40;
p(idx2) = 40*s(idx2).^0.3785;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
