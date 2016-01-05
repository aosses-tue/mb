function s = phon2sone(p)
% function s = phon2sone(p)
%
% 1. Description:
%
% 2. Stand-alone example:
%       p = [20 45];
%       s = phon2sone(p);
% 
% 3. Additional info:
%       See also: sone2phon
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 04/01/2016
% Last update on: 04/01/2016 
% Last use on   : 04/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx1 = find(p >= 40);
idx2 = find(p < 40);

s(idx1) = 2.^((p(idx1)-40)/10);
s(idx2) = (p(idx2)/40).^2.642;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
