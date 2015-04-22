% Combined outer and middle ear filter Filter for use with the DRNL
% Use this script after having generated the filters
% Morten Løve Jepsen, 2006

function y = OuterMiddleFilter(x);

load OuterFIR_44100CI.mat
load MiddleFIR_44100CI.mat

y_tmp = filter(outer_ear_fir_coeff,1,x);
y = filter(mid_ear_fir_coeff,1,y_tmp);

%eof
