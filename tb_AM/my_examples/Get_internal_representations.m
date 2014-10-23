function [outsig, fc, t] = Get_internal_representations(insig,fs,model)
% function [outsig, fc, t] = Get_internal_representations(insig,fs,model)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 22/10/2014
% Last update on: 22/10/2014 % Update this date manually
% Last use on   : 22/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    model = 'dau1996';
end

if nargin < 2
    error('Specify the sample rate of insig');
end

if strcmp(model,'dau1996a') % No overshoot limit
    [outsig fc] = dau1996apreproc(insig,fs);
end

if strcmp(model,'dau1996') % Overshoot limit
    [outsig fc] = dau1996preproc(insig,fs);
end

t = ( 1:size(outsig,1) )/fs;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
