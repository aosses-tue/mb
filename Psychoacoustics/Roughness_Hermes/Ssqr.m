function y=Ssqr(x)
% function y=Ssqr(x)
%
% 1. Description:
%       Cumulative sum of x-squared
% 
% 2. Example:
%       x = [1 2 3];  
%       y = Ssqr(x) % y = 1.*1+2.*2+3.*3 = 15
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 0; 
for a =1:1:length(x)
   z = z + (x(a)^2);
end

y =	z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end