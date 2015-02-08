% Received on 20150206 from AK (e-mail)
% 
% My description:
%       Quantisation per segment (50-ms length): looks for No of bits producing 
%       an SNR of approx. 30 dB

%% Quantize a file in blocks automatically

inputfile = 'trumpet.wav';      % Input file
outputfile = 'out.wav';         % Output file

Tb      = 50e-3;                % Block length (sec)
bmin    = 0;                    % Minimum number of bits
bmax    = 15;                   % Max number of bits
SNR     = 30;                   % Desired signal-to-noise ratio (dB)

% Load the file (fs = sampling freq)
[x,fs] = wavread(inputfile);

% Take the left channel to make sure the signal is mono
x=x(:,1);

% Number of samples
N = length(x);

% Number of samples per block
Ns = round(Tb*fs);

% Number of blocks
K = floor(N/Ns);

% Last block size
Nlast = mod(N,Ns);
if Nlast
    K=K+1;
else
    Nlast=Ns;
end

% Reserve memory for result
y = zeros(N,1);

% Reserve memory for bits per block
B = zeros(K,1);

% RMS values per block
R = zeros(K,1);

% Loop through the blocks
for k = 1:K
    
    % Last block:
    if k==K
        Nblock = Nlast;
    else
        Nblock = Ns;
    end
    
    % First sample of block
    N1 = (k-1)*Ns+1;
    
    % Last sample of block
    N2 = N1+Nblock-1;
    
    % Get data
    xb = x(N1:N2);
    
    % Calculate RMS (dB)
    R(k) = 10*log10(mean(xb.^2)+eps);

    S=0;
    b=bmin-1;
     
     while (S<SNR) && (b<bmax)
        % Increase bitrate
        b=b+1;
        % Quantize
        Q = 2^(b-1);
        yb = round(xb*Q)/Q;
        % Noise signal (=difference between original and coded)
        w = xb-yb;
        % RMS of the difference
        r = 10*log10(mean(w.^2));
        % Signal-to-noise
        S=R(k)-r;
     end
     
     % Store
    y(N1:N2) = yb;
    B(k)=b;
end

% Difference signal
d = y-x;

% Plot
figure
N = size(x,1);      % Number of samples
t = (0:N-1)'/fs;    % Time axis

tblock = buffer(t,Tb*fs,0);
tblock = tblock(1,:);

% Original
subplot(4,1,1);
plot(t,x);
axis([t(1) t(end) -1.1 1.1]);
xlabel('time (s)');
ylabel('Original');
title(inputfile)

% Coded
subplot(4,1,2);
plot(t,y);
axis([t(1) t(end) -1.1 1.1]);
xlabel('time (s)');
ylabel('Coded');

% Difference
subplot(4,1,3);
plot(t,d);
axis([t(1) t(end) min(d) max(d)]);
xlabel('time (s)');
ylabel('Difference');

% Mean squared error
e = sqrt(mean(d(:).^2));
fprintf('Root mean squared error: %.2e (%.2f dB)\n', e, 20*log10(e));

% Mean number of bits
fprintf('Mean number of bits: %.2f\n', mean(B));

subplot(4,1,4);
plot(tblock,B);

xlim([t(1) t(end)])
xlabel('Block number')
ylabel('Num. of bits')
disp('')