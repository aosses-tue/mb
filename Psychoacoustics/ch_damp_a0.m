function fgp_d = ch_damp_a0(fgrp, HL_ihc)
% function fgp_d = ch_damp_a0(fgrp, HL_ihc)
%
% 1. Description:
%       Attenuation due to outer & middle ear and inner hair cell hearing loss
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000 (new version with comments and examples on 06/01/2007)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 23/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0 = [ 0 0 0 0 0 0 0 0 0 0 -.2 -.5 -1.2 -2.1 -3.2 -4.6 -5.5 -5.6 -4.3 -2.5 -0.1 2.8 6.4 20.0];
fgp_d=zeros(length(fgrp(:,1)),24);  
for i = 1:length(fgrp(:,1))
   fgp_d(i,:)=fgrp(i,:)-a0-HL_ihc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
