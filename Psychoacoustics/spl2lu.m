function [LU1 SPL1] = spl2lu(SPL1)
% function [LU1 SPL1] = spl2lu(SPL1)
%
% 1. Description:
%       Converts SPL into loudness units LU, according to Fletcher1933
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/09/2014
% Last update on: 05/09/2014 % Update this date manually
% Last use on   : 05/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spl = -10:10:120;
lu = [0.015 1 13.9 97.5 360 975 2200 4350 7950 17100 38000 88000 215000 556000];

if nargin == 0
    SPL1 = spl;
end

LU1 = interp1(spl,lu,SPL1);

if nargout == 0
    figure;
    plot( SPL1,log(LU1) ,'o-');
    grid on
    xlabel('SPL [dB]')
    ylabel('Loudness Units in logarithm ln(LU)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
