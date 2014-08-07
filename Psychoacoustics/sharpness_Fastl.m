function [sharp] = sharpness_Fastl(loudspec)
% function [sharp] = sharpness_Fastl(loudspec)
% 
% Method FASTL (1991)
% Expression for weighting function obtained by fitting an equation to data 
% given in 'Psychoacoustics: Facts and Models' using MATLAB basic fitting 
% function
% sharp = sharpness [acum]
% 
% 1.1 EXAMPLE for a stationary sound (1 kHz sine-tone):
%       f       = 1000; % [Hz]
%       dur     = 1; % [s]
%       fs      = 44100;
%       [Sig, t]  = Create_sin(f,dur,fs);
%       Sig     = setdbspl(Sig,70);
%       [N_tot, loudspec, BarkAxis, LN] = Loudness_ISO532B_from_sound(Sig,fs);
%       sharp   = sharpness_Fastl(loudspec);
%       sprintf('sine RMS = %.1f [dB], estimated sharpness = %.2f [acum]', dbspl(Sig), sharp)
% 
% 1.2 EXAMPLE for a stationary sound (1kHz +2kHz sine-tone) lying in 
%             different critical bands:
%       f1      = 1000; % [Hz]
%       f2      = 2000;
%       dur     = 1; % [s]
%       fs      = 44100;
%       [Sig1, t]  = Create_sin(f1,dur,fs);
%       [Sig2, t]  = Create_sin(f2,dur,fs);
%       Sig     = setdbspl(Sig1+Sig2,70);
%       [N_tot, loudspec, BarkAxis, LN] = Loudness_ISO532B_from_sound(Sig,fs);
%       sharp   = sharpness_Fastl(loudspec);
%       sprintf('2-component sine RMS = %.1f [dB], estimated sharpness = %.2f [acum]', dbspl(Sig), sharp)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Claire Churchill
% Created in    : Sept 2004
% Downloaded on : 22/07/2014 (http://www.salford.ac.uk/computing-science-engineering/research/acoustics/psychoacoustics/sound-quality-making-products-sound-better/sound-quality-testing/matlab-codes)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 22/07/2014 % Update this date manually
% Last use on   : 22/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(loudspec);

gz(1:140) = 1;
z = 141:n;
gz(z) = 0.00012*(z/10).^4-0.0056*(z/10).^3+0.1*(z/10).^2-0.81*(z/10)+3.5;

z = 0.1:0.1:(n/10);

sharp = 0.11 * sum(loudspec.*gz.*z.*0.1) / sum(loudspec.*0.1);
% sharp = 0.11 * sum(loudspec.*gz.*z.*0.1) / sum(loudspec.*0.1 +eps); % add eps if necessary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
