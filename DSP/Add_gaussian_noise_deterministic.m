function [y noise] = Add_gaussian_noise_deterministic(x,mu,sigma)
% function [y noise] = Add_gaussian_noise_deterministic(x,mu,sigma)
%
% 1. Description:
%       It adds a deterministic gaussian noise with mean mu and standard 
%       deviation sigma to the input signal x. The array x should be a 
%       column vector.
%        
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also: Add_gaussian_noise.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 04/08/2015
% Last update on: 04/08/2015 % Update this date manually
% Last use on   : 04/08/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1);
M = size(x,2);

noise = normrnd(mu,sigma,N,1);
noise = repmat(noise,1,M);

y = x + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
