% Received on 20150206 from AK (e-mail)

%% Quantize a complete file

inputfile = 'trumpet.wav';      % Input file
outputfile = 'out.wav';         % Output file
Nbits = 10;                     % Number of bits

% Load the file (fs = sampling freq)
[x,fs] = wavread(inputfile);

% Take the left channel to make sure the signal is mono
x=x(:,1);

% Quantize
G = 2^(Nbits-1);
y = round(x*G)/G;

% Store
wavwrite(y,fs,16,outputfile);

% Play
sound(y,fs);

% Difference signal
d = y-x;

% Plot
figure
N = size(x,1);      % Number of samples
t = (0:N-1)'/fs;    % Time axis

% Original
subplot(3,1,1);
plot(t,x);
axis([t(1) t(end) -1.1 1.1]);
xlabel('time (s)');
ylabel('Original');

% Coded
subplot(3,1,2);
plot(t,y);
axis([t(1) t(end) -1.1 1.1]);
xlabel('time (s)');
ylabel('Coded');

% Difference
subplot(3,1,3);
plot(t,d);
axis([t(1) t(end) min(d) max(d)]);
xlabel('time (s)');
ylabel('Difference');

% Mean squared error
e = sqrt(mean(d(:).^2));
fprintf('Root mean squared error: %.2e (%.2f dB)\n', e, 20*log10(e));
