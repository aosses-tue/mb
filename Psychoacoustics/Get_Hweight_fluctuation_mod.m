function Hweight = Get_Hweight_fluctuation_mod(N,Fs)
% function Hweight = Get_Hweight_fluctuation_mod(N,Fs)
%
% 1. Description:
%       Only HPF
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/11/2014
% Last update on: 25/11/2014 % Update this date manually
% Last use on   : 25/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hweight     = zeros(1,N);
Hweight_bef	= zeros(1,N);

DCbins	= 0;

H2 = [	0       0
        0.2441  0
        0.7324	0.45
        1.7090	0.82
        1.9531  0.9
        3.1738	1
        12.0    1
        15.     0.7% 0.55
        16      0 
        30      0
        150     0]; %0];     % obtained from Sontacchi1999, Abbildung 5.4, understanding
                        % that 'BINS' are converted to frequency in such a way that
                        % bin 16th is approximately 16 Hz
                        
fH2 = H2(:,1); 

% last	= floor((16/Fs)*N) ;
last	= min( floor((150/Fs)*N) ,N/2);
k       = DCbins+1:1:last;
f       = (k-1)*Fs/N;

idx     = 1:length(f);

% Hhann = zeros(size(Hweight));
% Hhann(1:8192) = Hanning_half(8192);

% opts.fs = Fs;
% [yhan ydB fhan] = freqfft(Hhann',N/2,opts);

% fref = 4;
% y4Hz = interp1(fhan,abs(yhan),fref);

% idx = 1:length(fhan);

Hweight_bef(1,k) = interp1(fH2, H2(:,2),f(k - DCbins));

Hweight(1,k) = Hweight_bef(1,k);
% Hweight(idx) = Hweight(idx).*abs(yhan');

% y4Hz = interp1(fhan,Hweight(1:length(fhan)),fref);

% if ~isnan(y4Hz) & y4Hz ~= 0
%     Hweight = Hweight/y4Hz;
% end

if nargout == 0
    figure;
    plot( f(idx),Hweight(idx) );
    legend('Hweight')
    xlim([0 30])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
