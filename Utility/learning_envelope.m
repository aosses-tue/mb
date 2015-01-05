function y = learning_envelope(x)
% function y = learning_envelope(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 03/01/2015
% Last update on: 03/01/2015 % Update this date manually
% Last use on   : 03/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPart1 = 0;
bPart2 = 0; % Overlap of windows given a size and overlap
bPart3 = 1;

if bPart1
%%
f       = 1000;
fs      = 44100;
N           = 512; % samples of the time series
WinLength   = 1024; % samples of the window

dur = (N/fs)
[y, t]  = Create_sin(f,dur,fs);
plot(t,y); xlabel('time [s]'); ylabel('Amplitude')
    
w = Get_window('hamming',WinLength-1);
w(WinLength) = 0;

y = Zero_padding(y,(WinLength-N)/fs,fs);

freqz(y.*w,1,WinLength)
end

if bPart2
%%
[x fs] = wavread('stft_audiofile.wav');
stft(x,fs);
[xx xx xx opts] = stft(x,fs);
w = Get_window('hamming',opts.wlen);

total_samples = 2*fs;

cols = 10;
y = zeros(total_samples,cols);

opts.wlen = 4096;
opts.hope = 4096/2;

for i = 1:cols
    y(:,i) = Do_zeropadding_in_samples(w, opts.hope*(i-1), total_samples);
end

yf = transpose( sum(transpose(y)) );

figure; plot(y,'b'), hold on; plot(yf,'r');
end

if bPart3
%%

f           = 1000;
fs          = 44100;
N           = 8192; % samples of the time series
WinLength   = 8192; % samples of the window
fmod        = 550; % modulation frequency

dur = (N/fs)
[x, t]  = Create_sin(f,dur,fs);
    
w = Get_window('blackman',WinLength);

m           = 100-11;
option      = 'm';
start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
[y,env]     = ch_am(x,fmod,m,option,fs,start_phase);
yw          = y.*w;

etmp        = fft(yw)

close all
figure;
plot(t,yw); xlabel('time [s]'); ylabel('Amplitude')
title('Windowed signal to which envelope is going to be extracted')

% ei          = N*real(ifft(etmp));
% etmp_td     = abs(ei);
etmp_td     = abs(yw);
h0          = mean(etmp_td);

Hw_fs       = Get_Hweight_fluctuation_mod(N,fs);
Hw_fs       = transpose(Hw_fs);

Hw_r        = Get_Hweight_roughness(N,fs);
Hw_r        = transpose( Hw_r(17,:) );

Fei         = fft( etmp_td-h0 ); % changes the phase but not the amplitude
hBPi        = 2*real(  ifft( Fei.*Hw_r )  );
% hBPi        = 2*real(  ifft( Fei.*Hw_fs)  );

figure;
plot(t,etmp_td, t,hBPi);

xlabel('time [s]'); ylabel('Amplitude')
title('Envelope e_i')

hBPrms      = dw_rms(hBPi);
mdepht      = hBPrms / h0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
