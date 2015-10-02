function [detect,prob] = caspmdecide(mu,in_var,rule,numint)
% function [detect,prob] = caspmdecide(mu,in_var,rule,numint)
%
% Y = 1, signal is detected
%   
%   rule = 1 (x-down, 1 up procedures)
%   numint = 2, 3 or 4 for 2, 3 or 4-AFC
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/caspmdecide.php
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% 1. Stand-alone example:
%   mu = 45
%   in_var = 1.1;
%   rule = [1 2]; % 2-down, 1-up
%   numint = 3;
%   [detect prob] = caspmdecide(mu,in_var,rule,numint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    numint = 3; % 3-AFC
end
if nargin < 3
   rule = [1 2]; % 2-down, 1-up
end

if in_var <= 0
    error('CASP:mdecide', 'in_var must be > 0');
end

if rule(1) > 1
    error('CASP:mdecide','Only x-down, 1-up procedure have been implemented')
end
 
switch numint
    case 2 % 2-AFC
        prob = 1 - (erfc((((mu/in_var)*0.707) - 0    ) * sqrt(2)/2) / 2); % complementary error function
    case 3 % 3-AFC
        prob = 1 - (erfc((((mu/in_var)*0.765) - 0.423) * sqrt(2)/2) / 2);
    case 4 % 4-AFC
        prob = 1 - (erfc((((mu/in_var)*0.810) - 0.668) * sqrt(2)/2) / 2);
    otherwise
        error('CASP:mdecide', 'Only 2-, 3- and 4-AFC procedures are implemented');
end;

if rule(1) == 1 && max(prob) > (1 / (2 .^ (1/rule(2))))
    detect = 1; % signal heard
else
    detect = 0; % no signal heard
end;

%OLDFORMAT
