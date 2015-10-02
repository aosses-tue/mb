function Pcorrect = ideal_observer(SNRenv,parameters)
% function Pcorrect = ideal_observer(SNRenv_lin,parameters)
%
% 1. Description:
%       Copied from inline function in the script joergensen2013.
% 
% 2. Stand-alone example:
%       SNRenv = 2;
%       % parameters = [0.61 0.5 8000 0.6]; % Run the model with parameters for the CLUE material from Joergensen et al., (2013).
%       parameters = [0.61 0.5 3 0.6]; 
%       Pcorrect = ideal_observer(SNRenv,parameters); 
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/10/2015
% Last update on: 01/10/2015 
% Last use on   : 01/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin < 2
%     error('You have to specify the k,q,m,sigma_s parameters for the IdealObserver')
% end

k = parameters(1);
q = parameters(2);
m = parameters(3);
sigma_s = parameters(4);

% ---------- Converting from SNRenv to d_prime  --------------
d_prime = k*(SNRenv).^q;

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
Un = 1*norminv(1-(1/m));
m_n = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
Pcorrect = normcdf(d_prime,m_n,sqrt(sigma_s.^2+sig_n.^2))*100;

% Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
% In Signal Detection and Recognition by Human Observers,
% edited by John A. Swets (John Wiley & Sons, New York)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
