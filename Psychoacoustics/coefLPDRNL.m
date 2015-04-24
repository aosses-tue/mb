function [b,a]=coefLPDRNL(fc,fs);
% function [b,a]=coefLPDRNL(fc,fs);
% 
% 1. Description:
%       Designs a 2nd order low-pass filter
% 
% rev Morten Loeve Jepsen, 2.nov 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = pi*fc/fs;

C = 1/(1+sqrt(2)*cot(theta)+(cot(theta))^2);
D = 2*C*(1-(cot(theta))^2);
E = C*(1-sqrt(2)*cot(theta)+(cot(theta))^2);

b = [C, 2*C, C];
a = [1, D, E];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
