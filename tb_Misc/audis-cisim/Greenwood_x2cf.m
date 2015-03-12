function cf = Greenwood_x2cf(x);
%Calculates the center frequencies corresponding to a place along the Basilar Membrane
%Based upon Greenwoods formula
%
% cf = Greenwood_x2cf(x);
%
% x is the distance from the apical end of the Basilar Membrane expressd in mm
% cf is the resulting frequency in Hertz
%
% Copyright Tom Francart, 2009



% define constants of Greenwoods formula
A=165.4;
a=0.06;
k=1;

cf = A * (10.^(a*x) - k);
return;