function y = r20141121_DSP_and_FS(x)
% function y = r20141121_DSP_and_FS(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/11/2014
% Last update on: 20/11/2014 % Update this date manually
% Last use on   : 20/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 8192;
% w = Get_window('hanning',N);
w = hanning(N,'symmetric'); % I HAVE to use a symmetric window
w(N/2+1:end) = 0;

w = [w; zeros(N,1)];

[w2 corrdB] = Hanning_half(N);
K = N/2;

opts.fs = 44100;
[y  ydB  f] = freqfft(From_dB(corrdB)*w ,K,opts);
[y2 ydB2 f] = freqfft(From_dB(corrdB)*w2,K,opts);

t = ( 1:length(w) )/opts.fs;

figure; 
subplot(2,1,1)
plot(f,ydB); grid on, hold on
plot(f,ydB2,'r');
legend('half hanning','half hanning, script')
xlim([0 30])
grid on

subplot(2,1,2)
plot(f,abs(y)); grid on, hold on
plot(f,abs(y2),'r');
legend('half hanning','half hanning, script')
xlim([0 30])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
