function y = r20151204_Antoine_lecture
% function y = r20151204_Antoine_lecture
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 04/12/2015
% Last update on: 04/12/2015 
% Last use on   : 04/12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
dir = Get_TUe_data_paths('piano');

bDoPart1 = 0;
bDoPart2 = 1;

file = [dir '01-Chabassier\SONS\Cd5\pressionexpe.wav'];
titlelabel = 'Cd5';

if bDoPart1
    
    [x fs] = Wavread(file);
    t = ( 1:length(x) )/fs;

    % figure;
    % plot(insig);

    Ni = 3191; % max of the waveform (manually computed)
    Nf = 18080; % end of analysis period

    insig = x(Ni:Nf);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slide 11/
    N = 8192;
    K = N/2;
    windowtype = 'hanning';

    figure;
    subplot(2,2,1)
    freqfft2(insig,K,fs,windowtype);     
    xlim([0 10000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))

    subplot(2,2,3)
    freqfft2(insig,K,fs,windowtype);     
    xlim([3000 5000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))

    %%%
    N = 65536;
    K = N/2;

    subplot(2,2,2)
    freqfft2(insig,K,fs,windowtype);     
    xlim([0 10000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))

    subplot(2,2,4)
    freqfft2(insig,K,fs,windowtype);     
    xlim([3000 5000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slide 14/
    file = [dir '01-Chabassier\SONS\F3\pressionsimu.wav'];
    titlelabel = 'F3';

    [x fs] = Wavread(file);

    % figure;
    % plot(x);

    Ni = 4380; % max of the waveform (manually computed)
    Nf = 46000; % end of analysis period

    insig = x(Ni:Nf); 

    N = 16384;
    K = N/2;

    windowtype   = 'rectangular';
    [y0, y0dB f] = freqfft2(insig,K,fs,windowtype);     

    windowtype   = 'hanning';
    [y1, y1dB]   = freqfft2(insig,K,fs,windowtype);     

    windowtype   = 'flattop';
    [y2, y2dB]   = freqfft2(insig,K,fs,windowtype);     

    figure;
    plot(f,y0dB, f,y1dB, f,y2dB)
    xlim([0 2000])
    ylim([0 80])
    legend('rect','hann','flattop');
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz], use of different windows',titlelabel,N,fs))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slide 15/
    A = 1;
    f1 = 100;
    f2 = 200.5;
    fs = 1000;
    Nsam = 1000;
    x1 = Create_sin(f1,Nsam/fs,fs);
    x2 = Create_sin(f2,Nsam/fs,fs);
    xtot = x1+x2;

    windowtype   = 'rectangular';

    N = 1024; % resolution of 1 Hz at fs = 1000 Hz
    K = N/2;
    [y0, y0dB, f0] = freqfft2(xtot,K,fs,windowtype);     

    N = 2048; % resolution of 1 Hz at fs = 1000 Hz
    K = N/2;
    [y1, y1dB, f1]   = freqfft2([xtot; zeros(1024,1)],K,fs,windowtype);     

    figure;
    plot(f0,abs(y0),f1,abs(y1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slide 18/

    file = [dir '01-Chabassier\SONS\Cd5\pressionexpe.wav'];
    titlelabel = 'Cd5';

    [x fs] = Wavread(file);
    t = ( 1:length(x) )/fs;

    % figure;
    % plot(insig);

    Ni = 3191; % max of the waveform (manually computed)
    Nf = 18080; % end of analysis period

    insig = x(Ni:Nf);

    N = 65536; % N-point FFT
    K = N/2;
    windowtype = 'hanning';

    figure;
    subplot(2,1,1)
    freqfft2(insig,K,fs,windowtype);     
    xlim([3000 5000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))

    %%%
    N = 16384;
    K = N/2;
    [STFT f t] = stft(insig,fs,N,K,75,'hanning');

    subplot(2,1,2)
    plot(f,20*log10(abs(STFT(:,1)))); grid on
    xlim([3000 5000])
    title(sprintf('%s - N=%.0f, fs=%.0f [Hz]',titlelabel,N,fs))

end

if bDoPart2
    [x fs] = Wavread(file);
    t = ( 1:length(x) )/fs;

    % figure;
    % plot(insig);

    % L = 50;
    N = 300; % N has to be greater than 5*L
    
    Ni = 3191; % max of the waveform (manually computed)
    Nf = Ni+N-1; % 18080; % end of analysis period

    insig = x(Ni:Nf);
    
    p = 180;
    X0 = il_get_X0(insig,p,'X0');
    X1 = il_get_X0(insig,p,'X1');
    X0plus = pinv(X0);
    Matr = X0plus*X1;
    
    lambda = eig(Matr);
    
    idxposim = find(imag(lambda)>0);
    
    L = length(idxposim);
    
    %%% Frequency and damping factors:
    fi = ( atan2(imag(lambda(idxposim)),real(lambda(idxposim))) )/(2*pi);
    Fi = fi*fs;
    alphai = -log10(abs(lambda(idxposim)));
    Alphai = alphai*fs;
    
    M = N;
    %%% Form matrix T
    T = nan(M,2*L);
    n = 1;
    
    for n = 1:M
        nu_i(n,:)   = transpose(exp(-n*alphai) .* cos(2*j*pi*n*fi));
        beta_i(n,:) = transpose(exp(-n*alphai) .* sin(2*j*pi*n*fi));
    end
    
    idxs = 1:2:2*L;
    T(1,idxs) = 1;
    
    idxs = 2:2:2*L;
    T(1,idxs) = 0;
    
    for n = 1:L
        if mod(n,2) == 1
            T(2:end,2*n-1) = nu_i(2:end,n);
        else
            T(2:end,2*n  ) = beta_i(2:end,n);
        end
    end
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

function A = il_get_X0(insig,p,type)

N = length(insig);

sizeN = N-p;
sizeM = p;

A = nan(sizeN,sizeM); % memory allocation

switch type
    case 'X0'
        increi = 1;
    case 'X1'
        increi = 2;
end

for i = 1:sizeM
    
    idx1 = p-i+increi;
    idx2 = N-i+increi-1;
    A(:,i) = insig(idx1:idx2);
    
end
disp('')