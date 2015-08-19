function out = Get_template(ir1,ir2,fs)
% function out = Get_template(ir1,ir2,fs)
%
% 1. Description:
%       ir1 - internal representation of signal 1
%       ir2 - internal representation of signal 2
% 
%           ir1 - related to 'Noise alone'
%           ir2 - related to 'Suprathreshold signal' (signal well above threshold)
%           fs - Sampling frequency: relevant for the normalisation process.
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 23/10/2014
% Last update on: 30/04/2015 
% Last use on   : 19/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
	fs = 44100;
end

idx = 1:size(ir1,2);

[xx xx out] = Normalise_signal(ir2(:,idx)-ir1(:,idx),fs);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
