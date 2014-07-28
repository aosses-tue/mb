function f = Get_OB_freqs(BandsPerOctave)
% function f = Get_OB_freqs(BandsPerOctave)
%
% 1. Description:
%   
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       BandsPerOctave = 1; % Octave bands
%       f = Get_OB_freqs(BandsPerOctave);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 3/7/2014
% Last update on: 3/7/2014 % Update this date manually
% Last used on  : 3/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    BandsPerOctave = 1;
end

fs  = 44100;    % Sampling frequency
N   = 8;        % Filter Order
F0  = 1000;     % Center Frequency (Hz)
filtSpecs   = fdesign.octave(BandsPerOctave,'Class 1','N,F0',N,F0,fs);
f   = validfrequencies(filtSpecs);

% for i=1:length(f)
    % filtSpecs.F0 = f(i);
    % Hd3(i) = design(filtSpecs,'butter');
    % Get_filter_specs(Hd3(i));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end