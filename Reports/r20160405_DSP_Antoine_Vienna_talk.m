function y = r20160405_DSP_Antoine_Vienna_talk
% function y = r20160405_DSP_Antoine_Vienna_talk
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 05/04/2016
% Last update on: 05/04/2016 
% Last use on   : 05/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
[x fs] = Wavread(file);

t = ( 1:length(x) )/fs;

x_pa = 2*x;
pr = abs(x_pa);
pr = Apply_IIR_Butter(pr,fs,20,'low',4);

bDoAlles = 0;

if bDoAlles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1:
figure;
plot(t,To_dB(abs(x_pa)/2e-5)); hold on

plot(t,To_dB(2*pr/2e-5),'r','LineWidth',2);
ylim([35 110])
xlabel('Time [s]')
ylabel('Sound Pressure Level [dB]')
title('Envelope Cd5')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2:

Nstart = 3190;
Nend   = Nstart + round(0.5*fs)-1;

x_sec = x_pa(Nstart:Nend);

N = 8192;
K = N/2;
[y ydB f] = freqfft2(x_sec,K,fs,'rectangular');

figure;
plot(f,ydB);

xlim([0 10000])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
title(sprintf('Pr1 - Cd5 - N=%.0f - fs = %.0f [kHz]',N,fs/1000))

figure;
plot(f,ydB);

xlim([3000 5000])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
title(sprintf('Pr1 - Cd5 - N=%.0f - fs = %.0f [kHz]',N,fs/1000))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 2^16;
K = N/2;
[y ydB f] = freqfft2(x_sec,K,fs,'rectangular');

figure;
plot(f,ydB);

xlim([0 10000])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
title(sprintf('Pr1 - Cd5 - N=%.0f - fs = %.0f [kHz]',N,fs/1000))

figure;
plot(f,ydB);

xlim([3000 5000])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
title(sprintf('Pr1 - Cd5 - N=%.0f - fs = %.0f [kHz]',N,fs/1000))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df = 1; % Hz
N = round(fs/df);
K = N/2;
[y ydB f] = freqfft2(x_sec,K,fs,'rectangular');

[yhan ydBhan] = freqfft2(x_sec,K,fs,'hanning');

figure; 
plot(f,ydB); hold on;
plot(f,ydBhan,'r')
legend('rect','hanning')
xlim([400 4200])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Zero padding

f1 = 100;
f2 = 200.5;
A1 = 1;
A2 = 1;
fs = 1000;
N = 1000;
y1 = A1*Create_sin(f1,N/fs,fs);
y2 = A2*Create_sin(f2,N/fs,fs);
yt = y1 + y2;

[o odB f] = freqfft2(yt,N/2,fs,'rectangular');
figure;
subplot(1,2,1)
plot(f,abs(o));
 
[o2 odB2 f] = freqfft2([yt; zeros(N,1)],N,fs,'rectangular');
%[o2 odB2 f] = freqfft2(yt,N,fs,'rectangular');
subplot(1,2,2)
plot(f,abs(o2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 4000;
A = [1 0.7 1.2 0.4];
f = [100 131 297 525];
Tau = [0.1 0.047 0.022 0.01]; % s
phi = [0 0 0 0]; 

sig     = 1./Tau; % damping factor in 1/s
alpha   = sig/fs;
fd = f/fs;

n       = transpose(1:fs/2);
insig   = zeros(size(n));

for i = 1:length(A)
    insig = insig + A(i)*exp(-alpha(i)*n).*cos(2*pi*n*fd(i)+phi(i));
end

N   = 64;
df  = 62.5;
fs  = round(df*N);

K   = N/2;

[y ydB f] = freqfft2(insig,K,fs,'rectangular');
figure;
plot(f,abs(y)/max(abs(y)));
xlabel('Frequency [Hz]')
ylabel('Magnitude [normalised]')
title(sprintf('%.f sinus - FFT - N = %.0f',length(A),N));

L = 4;
p = 2*L;
[outsig fi ai alpha_i Phi_i L] = Get_ESPRIT_analysis(insig,p,N,fs);

Tau_est = 1./alpha_i;

figure; 
plot(insig+5); hold on
xlabel('Time [samples]')
ylabel('Amplitude')

plot(outsig,'r')
legend('original','synthetic')
title(sprintf('%.f sinus - ESPRIT reconstruction - N = %.0f',L,N));

end
