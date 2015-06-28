function Figures_Psychoacoustics_tex(fig_label)
% function Figures_Psychoacoustics_tex(fig_label)
%
% 1. Description:
%
% 2. Stand-alone example:
%       fig_label = 'fig_a0.b';
%       Figures_Psychoacoustics_tex(fig_label);
% 
%       fig_label = 'fig_a0.c';
%       Figures_Psychoacoustics_tex(fig_label);
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 26/06/2015
% Last update on: 26/06/2015 % Update this date manually
% Last use on   : 26/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    fig_label = 'fig_a0.b';
end

switch fig_label
    case 'fig_a0.b'
        a0tab =	Get_psyparams('a0tab');
        z = a0tab(:,1); 
        a0values = a0tab(:,2);
        f = bark2hz(z); 
        f(1) = 0.5; 

        HTres = Get_psyparams('HTres'); 
        zh = HTres(:,1); 
        HTresvalues = HTres(:,2);fh = bark2hz(zh); 

        figure;
        semilogx(f , a0values); grid on, hold on
        semilogx(fh,-HTresvalues,'r');
        axis([20 20000 -50 10])
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');

        legend('factor a0 from FF','-Hearing threshold')
        
    case 'fig_a0.c'
        % Outer and middle ear combined bandpass filter
        % (Pflueger, Hoeldrich, Riedler, Sep 1997)
        % Highpass component
        b = 0.109*[1 1];
        a = [1 -2.5359 3.9295 -4.7532 4.7251 -3.5548 2.139 -0.9879 0.2836];
        % Lowpass component
        d = [1 -2 1];
        c = [1 -2*0.95 0.95^2];
        num2 = conv(b, d);
        den2 = conv(a, c);
        N = 8192;
        K = N/2;
        fs = 44100;
        [h w f2] = freqz(num2,den2,K);
        a0values = 20*log10(abs(h));
        f = w/pi * fs/2;
        figure;
        semilogx(f , a0values); grid on, hold on
        axis([20 20000 -50 10])
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
