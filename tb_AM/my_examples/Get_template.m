function out = Get_template(ir1,ir2,setup)
% function out = Get_template(ir1,ir2,setup)
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
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 30/04/2015 % Update this date manually
% Last use on   : 05/08/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    setup = [];
end

setup = ef(setup,'idx',1:size(ir1,2));
idx = setup.idx;

disp('Make sure ir2 is N+Suprathreshold signal')
setup   = Ensure_field(setup,'fs',44100);
fs      = setup.fs;
out     = Normalise_signal(ir2(:,idx)-ir1(:,idx),fs);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
