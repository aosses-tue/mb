function y = Round(x,numDecimals)
% function y = Round(x,numDecimals)
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    numDecimals = 4;
end

y = round( x*1*10^(numDecimals))/10^(numDecimals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end