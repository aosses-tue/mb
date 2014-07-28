function ysmooth = smoothfft(y,f, info)
% function ysmooth = smoothfft(y,f, info)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 09/05/2014
% Last update: 21/05/2014 % Update this date manually
% Last used: 21/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 % info does not exist
    info = [];
end

info = Ensure_field(info,'typeplot',1);     % 1 = semilogx
                                            % 2 = linear scale

ysmooth = [];
K = length(f);

bTruncated = 0;
N = 51; % use an odd number

% y = [zeros((N-1)/2,size(y,2)); y];

for i = 1:K
    
%     dist    = ( 2^octaves*f(i) )/2; % back and forth
%     idx     = find(f < f(i)+dist & f >= f(i)-dist);
%     N       = length(idx);
    
    Ndist   = (N-1)/2;
    
    window  = zeros(K,1);
    
%     window_tmp = rectwin(N);
    
    liminf = i-Ndist;
    limsup = i+Ndist;
    
    if liminf < 1
        liminf = 1;
        bTruncated = 1;
    end
    if limsup > length(y)
        limsup = length(y);
        bTruncated = 1;
    end
    
    window(liminf:limsup) = 1;

    for j = 1:size(y,2)
        ysmooth(i,j) = sum(y(:,j).*window)/N;
        if bTruncated == 1
            ysmooth(i,j) = NaN;
        end
    end
    bTruncated = 0;
end

if nargout == 0
    
    info.fs = max(f)*2;
    figure;
    switch info.typeplot
        case 1
            semilogx(f,ysmooth), grid on
            xlim([20 info.fs/2])
            xlabel(['Frequency [Hz], (fs = ' num2str(info.fs) ' [Hz])'])
        case 2
            plot((1:K)/K,ysmooth), grid on
            xlabel('Normalised frequency (x\pi rad/sample)')
    end
    ylabel('Magnitude [dB]')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF - ' mfilename])