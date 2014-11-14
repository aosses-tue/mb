function Hweight = Get_Hweight_fluctuation(N,Fs)
% function Hweight = Get_Hweight_fluctuation(N,Fs)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/11/2014
% Last update on: 12/11/2014 % Update this date manually
% Last use on   : 12/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hweight	= zeros(1,N);

DCbins	= 0;

H2 = [	0       0
        0.2441  0
        0.7324	0.45
        1.7090	0.82
        1.9531  0.9
        3.1738	1
        12.9395 1
        14.4043 0.91
        15.3809	0.55
        15.8691 0.25  
        16      0];     % obtained from Sontacchi1999, Abbildung 5.4, understanding
                        % that 'BINS' are converted to frequency in such a way that
                        % bin 16th is approximately 16 Hz
                        
fH2 = H2(:,1); 

last	= floor((16/Fs)*N) ;
k       = DCbins+1:1:last;
f       = (k-1)*Fs/N;

Hweight(1,k) = interp1(fH2, H2(:,2),f(k - DCbins));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
