function ild = Get_ILD(x,y)
% function ild = Get_ILD(x,y)
%
% 1. Description:
%   ild = Get_ILD(x);   % x is an stereo file, Ch 1 = L; Ch 2 = R
%   ild = Get_ILD(x,y); % x, y are mono files, x = L; y = R
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example (function called inside lindemann1986.m):
%   demo_lindemann1986;
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 14/5/2014
% Last update: 14/5/2014 % Update this date manually
% Last used: 14/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    y = x(:,2);
    x = x(:,1);
end
ild = dbspl(y)-dbspl(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])