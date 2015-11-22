function outsig = Do_pitch_stretch(insig,fs,factor,mode)
% function y = Do_pitch_stretch(insig,fs,factor,mode)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 18/11/2015
% Last update on: 18/11/2015 
% Last use on   : 18/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    mode = 'semitone';
end

if nargin < 3
    factor = 5;
end

switch mode
    case 'semitone'
        num = round(100*(1/ (2^(factor/12))));
    case 'percentage'
        num = 100-round(factor);
end
        
den = 100;
timesfaster = num/den; % 0.75
outsig   = pvoc(insig,timesfaster,1024);

outsig = resample(outsig,num,den); % NB: 0.8 = 4/5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
