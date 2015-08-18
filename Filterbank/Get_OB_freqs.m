function f = Get_OB_freqs(BandsPerOctave,fmin,fmax)
% function f = Get_OB_freqs(BandsPerOctave,fmin,fmax)
%
% 1. Description:
%   
% 2. Stand-alone example:
%       BandsPerOctave = 1; % Octave bands
%       f = Get_OB_freqs(BandsPerOctave);
% 
%       % To find 1/3 octave bands between 250 and 1300 Hz
%       f = Get_OB_freqs(3,250,1300);
% 
% 3. Additional info:
%   Tested cross-platform: No 
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/07/2014
% Last update on: 17/08/2015 
% Last used on  : 17/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    BandsPerOctave = 1;
end

if BandsPerOctave == 1
    k = -5:4;
elseif BandsPerOctave == 3
    k = (-15:12)/3;
end

f = 1000*2.^(k);

if nargin >= 2
    idx = find(f<fmin);
    f(idx) = [];
end

if nargin >= 3
    idx = find(f>fmax);
    f(idx) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end