function x = Greenwood_cf2x(cf);
%Calculates the coordinates along the Basilar Membrane corresponding to a Center Frequency
%Based upon Greenwoods formula
%
% x = Greenwood_cf2x(cf);
%
% x is the distance from the apical end of the Basilar Membrane expressd in mm
% cf is the center frequency in Hertz
%
% See also Greenwood_x2cf
%
% Copyright Tom Francart, 2009



% define constants of Greenwoods formula
A=165.4;
a=0.06;
k=1;

x=log10(cf/A+k)/a;
return;