function r20160509_update_BIBD(n)
% function r20160509_update_BIBD(n)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 09/05/2016
% Last update on: 09/05/2016 
% Last use on   : 09/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 9; % Treatments
% n = 7, b = 35 x 20 s = 11.6 min per block, b_red = 14 (lambda = 2), need 2.5 persons
% n = 9, b = 84 x 20 s = 11.6 min per block, b_red = 24 (lambda = 2), need 3.5 persons

lambda = 1;

b_trials = n*(n-1)*(n-2)/6;

b_reduced_by = (n-2)/lambda;

b_trials_person = b_trials / b_reduced_by;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
