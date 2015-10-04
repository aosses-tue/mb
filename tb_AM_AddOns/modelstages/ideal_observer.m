function [detect Pcorrect] = ideal_observer(d_prime,sigma_s,m,rule)
% function [detect Pcorrect] = ideal_observer(d_prime,sigma_s,m,rule)
%
% 1. Description:
%       Copied from inline function in the script joergensen2013.
% 
% 2. Stand-alone example:
%       % m = 8000; parameters = [0.61 0.5 0.6]; % Run the model with parameters for the CLUE material from Joergensen et al., (2013).
%       m = 3; parameters = [0.61 0.5 0.6];
%       SNRenv = 2;
%       k = parameters(1);
%       q = parameters(2);
%       sigma_s = parameters(3);
%       d_prime = k*(SNRenv).^q;
%       [detect Pcorrect] = ideal_observer(d_prime,sigma_s,m); 
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/10/2015
% Last update on: 02/10/2015 
% Last use on   : 02/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    rule = 2; % 2-down, 1-up
end
    
if nargin < 3
    m = 3;
end

% Converting from d_prime to Percent correct, Green and Birdsall (1964)----
Un = 1*norminv(1-(1/m));
m_n = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
Pcorrect = normcdf(d_prime,m_n,sqrt(sigma_s.^2+sig_n.^2));

% Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
% In Signal Detection and Recognition by Human Observers,
% edited by John A. Swets (John Wiley & Sons, New York)

if max(Pcorrect) > (1 / (2 .^ (1/rule)))
    detect = 1; % signal heard
else
    detect = 0; % no signal heard
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
