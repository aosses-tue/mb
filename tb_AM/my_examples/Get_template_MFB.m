function out = Get_template_MFB(ir1,ir2,fs)
% function out = Get_template_MFB(ir1,ir2,fs)
%
% 1. Description:
%       ir1 - internal representation of signal 1
%       ir2 - internal representation of signal 2
% 
%           ir1 - related to 'Noise alone'
%           ir2 - related to 'Suprathreshold signal' (signal well above threshold)
%           setup.fs - Sampling frequency: relevant for the normalisation process.
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/08/2015
% Last update on: 11/08/2015 
% Last use on   : 11/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    fs = 44100;
    warning('Default sampling frequency being used')
end

[N M]   = size(ir1);
ir1 = ir1(:);
ir2 = ir2(:);

out     = Normalise_signal(ir2-ir1,fs);
out     = reshape(out,N,M);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
