function courseASP_lec2(x)
% function courseASP_lec2(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/01/2015
% Last update on: 03/01/2015 % Update this date manually
% Last use on   : 03/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 2T1.1

factor = From_dB(-48.8015);

if isunix
    main_folder = '~/Documenten/Documenten-TUe/09-Training+activities/20141001-audio-signal-processing/sms-tools/';
end

dir_sounds = [main_folder 'sounds' delim];

filename = [dir_sounds 'oboe-A4.wav'];

[x, fs] = Wavread(filename);
pin = 5000+1;
N = 512;
w = Get_window('hamming',N-1,'symmetric');
w(512) = 0;

% ws = Get_window('hanning',N,'symmetric');
% wp = Get_window('hanning',N,'periodic');
% figure; plot(ws), hold on; plot(wp,'r')

x1 = x(pin-N/2:pin+N/2-1);
t1 = -N/2:N/2-1;

figure(1)
subplot(3,1,1)
plot(t1,x1)
legend('Amplitude')

xw = x1.*w; % windowed input signal

X = fft(xw,N);
X1 = X(1:N/2);
fbin = 1:N/2;

phase = angle(X);
phaseu = unwrap(angle(X));
 
mX = 20*log10( abs(X1)*factor ); % mX = 20*log10( abs(X1)*factor );

subplot(3,1,2)
plot(fbin,mX)
ylabel('Amplitude [dB]')

subplot(3,1,3)
plot(fbin,(phase(1:N/2))/pi, ...
     fbin,(phaseu(1:N/2))/pi )
legend('normal','unwrapped')
ylabel('Phase [rad]')

% 48.8015 dB difference with Python

%% Figure 2T2
% complex-sinewaves.py

%% Figure 2T2.1

filename = [dir_sounds 'violin-B3.wav'];

[x, fs] = Wavread(filename);
pin = 5000+1;
N = 1024;
w = Get_window('hamming',N-1,'symmetric');
w(N) = 0;

% ws = Get_window('hanning',N,'symmetric');
% wp = Get_window('hanning',N,'periodic');
% figure; plot(ws), hold on; plot(wp,'r')

x1 = x(pin-N/2:pin+N/2-1);
t1 = -N/2:N/2-1;

figure(2)
subplot(3,1,1)
plot(t1,x1)
legend('Amplitude')

xw = x1.*w; % windowed input signal

X = fft(xw,N);
X1 = X(1:N/2);
fbin = 1:N/2;

phase = angle(X);
phaseu = unwrap(angle(X));
 
mX = 20*log10( abs(X1)*factor ); % mX = 20*log10( abs(X1)*factor );

subplot(3,1,2)
plot(fbin,mX)
ylabel('Amplitude [dB]')

subplot(3,1,3)
plot(fbin,(phase(1:N/2))/pi, ...
     fbin,(phaseu(1:N/2))/pi )
legend('normal','unwrapped')
ylabel('Phase [rad]')


% dft-complex-sine-1.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
