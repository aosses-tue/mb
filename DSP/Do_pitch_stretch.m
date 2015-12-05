function outsig = Do_pitch_stretch(insig,factor,mode)
% function y = Do_pitch_stretch(insig,factor,mode)
%
% 1. Description:
%       The input signal insig is increased/decreased in pitch according to
%       the parameters factor and mode. A positive factor will increase the 
%       pitch while a negative factor will decrease it. The variable mode
%       determines whether the increase will be done in semitones or in
%       percentage. For the first option mode has to be set to 'semitone' and
%       for the latter to 'percentage'.
%       The sampling frequency of the input signal is independent of this
%       relative increase/decrease in pitch (percentage).
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20151119_piano_sounds
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 18/11/2015
% Last update on: 18/11/2015 
% Last use on   : 02/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    mode = 'semitone';
end

if nargin < 2
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

outsig = resample(outsig,num,den); % if den > num outsig the pitch shift is going up (asuming that fs will be constant)
                                   % if den < num outsig the pitch shift is going down

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
