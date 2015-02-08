% Received on 20150206 from AK (e-mail)

%% Quantize a file in blocks automatically

inputfile   = 'trumpet.wav';    % Input file
outputfile  = 'out.wav';        % Output file

Tb  = [0   1  2  3  4];         % Block start times (sec)
B   = [12 10 12 10 10];         % Bits per block

% Load the file (fs = sampling freq)
[x,fs] = wavread(inputfile);

% Take the left channel to make sure the signal is mono
x=x(:,1);

% Number of samples
N = length(x);

% Number of blocks
K = length(Tb);

% Reserve memory for result
y = zeros(N,1);

% Block lengths
L = zeros(1,K);

% Loop through the blocks
for k = 1:K
    
    % First sample of block
    N1 = floor(Tb(k)*fs)+1;
    
    % Last block:
    if k==K
        N2 = N;
    else
        N2 = floor(Tb(k+1)*fs);
    end
    
    % Get data
    xb = x(N1:N2);
    
    % Quantize
    Q = 2^(B(k)-1);
    yb = round(xb*Q)/Q;
    
    % Store
    y(N1:N2) = yb;
    
    % Block length (sec)
    LN = N2-N1+1;
    L(k) = LN/fs;
end

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
axis([t(1) t(end) 1.1*min(d) 1.1*max(d)]);
xlabel('time (s)');
ylabel('Difference');

% Mean squared error
e = sqrt(mean(d(:).^2));
fprintf('Root mean squared error: %.2e (%.2f dB)\n', e, 20*log10(e));

% (Weighted) mean number of bits
T = N/fs;
MB = sum(L.*B)/T;
fprintf('Weighted mean number of bits: %.2f\n', MB);
