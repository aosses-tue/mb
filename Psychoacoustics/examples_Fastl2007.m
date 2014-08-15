function examples_fastl2007
% function examples_fastl2007
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% AM tone
fc = 2000;
T = 1/fc;
modfreq = 166.67*2;
modT = 1/modfreq;

filename = 'D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Fastl2007-psychoacoustics\Extracted\track_04_t03_AM.wav';
[x fs] = Wavread(filename);
t = (1:length(x))/fs;

xmod = Create_sin(modfreq,max(t),fs);

ti = 1;
figure;
plot(t,x,'b',t,xmod,'r')
ylabel('amplitude')
xlabel('t [s]')
xlim([ti ti+5*modT])

figure
opt.fs = fs;
freqfft(x,4096,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
