function [outsig f0] = Do_pitch_stretch(insig,factor,mode,fs)
% function [outsig f0] = Do_pitch_stretch(insig,factor,mode,fs)
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
% Last use on   : 22/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    mode = 'semitone';
end

if nargin < 2
    factor = 5;
end

den = 10000;
switch mode
    case 'semitone'
        num = round(100*(1/ (2^(factor/12))));
    case 'percentage'
        num = den-round(den*factor/100);
end
        
timesfaster = num/den; % 0.75
N       = 4096; % N = 1024 was giving some phase problems
outsigtmp  = pvoc(insig,timesfaster,N); % produces a time-stretched output signal (same pitch)

% num = round( factor*fs );
% den = fs;
outsig  = resample(outsigtmp,num,den); % if den > num outsig the pitch shift is going up (asuming that fs will be constant)
                                    % if den < num outsig the pitch shift is going down
if nargout == 0 | nargout == 2
    
    f0(1) = il_get_ESPRIT(insig,fs);
    f0(2) = il_get_ESPRIT(outsig,fs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f0 = il_get_ESPRIT(insig,fs,fmax)

if nargin < 3
    fmax = 2000;
end

N   = round(0.1*fs); % N has to be greater than 5*L, arbitrarily chosen
[insig_max Ni] = max(abs(insig)); % max of the waveform (manually computed)
Ni = Ni + round(0.020*fs); % 20-ms after the max
Nf  = Ni+N-1; 
insig_t = insig(Ni:Nf);
p = 180;
[out_ Fi Ai] = Get_ESPRIT_analysis(insig_t,p,N,fs);

idxx = find(Fi < fmax);
Fi = Fi(idxx); 
Ai = Ai(idxx);

[Amax idxmax] = max(Ai); % first 10 partials
idxFi = idxmax;
Fit = Fi;

factor2div = 1;

Ait = Ai(idxFi);
Fit = Fi(idxFi)/factor2div; % if factor2divide ~= 0, then 'virtual pitch' (when using ESPRIT) 
[Amax idxmax] = max(Ait); % first 10 partials
f0 = Fit(idxmax);

if factor2div == 1
    fprintf('\tf0 = %.3f [Hz]\n\n',f0);
else
    fprintf('\tf0 = %.3f [Hz] - ''virtual pitch'' \n\n',f0);
end