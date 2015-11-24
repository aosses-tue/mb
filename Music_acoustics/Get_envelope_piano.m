function outsig = Get_envelope_piano(insig,fs)
% function outsig = Get_envelope_piano(insig,fs)
%
% 1. Description:
%       Applies a 4th-order Butterworth filter centred at 20-Hz.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20151119_piano_sounds.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/11/2015
% Last update on: 24/11/2015 
% Last use on   : 24/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yin    = abs(hilbert( insig ));
[b, a] = butter(4,20/(fs/2),'low');

ylp    = filtfilt(b,a,yin); % linear phase implementation

% Aslow = 8; Apeak = 6; % Parameters used by RK
Aslow = 2; 
Apeak = 0; % Momentarily onset enhancement is by-passed!

ydiff = yin - Aslow*ylp;  % Makes Aslow*ylp well below 0 (except for the onset)
                        % figure; plot(ydiff)
ydiff = max(ydiff,0);
yout = ydiff*Apeak;

outsig = ylp; % yin + yout;

if nargout == 0
    % figure; freqz(b,a,4096);
    figure;
    subplot(2,1,1)
    plot(yin); hold on
    plot(yout,'r');
    ha = gca;
    % sound(ytot,fs)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
