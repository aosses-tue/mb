function [outsig Fi Ai sigma_i Phi_i L] = Get_ESPRIT_analysis(insig,p,M,fs)
% function [outsig Fi Ai sigma_i Phi_i L] = Get_ESPRIT_analysis(insig,p,M,fs)
%
% 1. Description:
%       outsig is a synthesised signal by approximating the insig signal
%       using the pencil matrix method.
% 
%       outsig = sum_{i=1}^L A_i \exp{-\alpha_i n} \cdot cos(2*pi*n*f_i + \phi_i)
% 
%       Fi = fi * fs;   % [Hz], fi is expressed in [decrement/sample]
%       sigma_i = alpha_i * fs; % [1/s] damping factor
% 
%       Input parameters:
%        - M samples of the input signal are used
%        - p is the pencil parameter
% 
%       Output parameters:
%        - outsig will be approximated as the sum of L damped sinuoids. With
%           L calculated from N and p. 2*L <= p <= N-2*L
%        - Fi [Hz], frequency of the sinusoids
%        - Ai, amplitude of the sinusoids
%        - Alpha_i [1/s], damping factor. Tau = 1./Alpha_i;
% 
% 2. Stand-alone example:
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [x fs] = Wavread(file);
%       N   = 2400; % N has to be greater than 5*L, arbitrarily chosen
%       Ni  = 3191; % max of the waveform (manually computed)
%       Nf  = Ni+N-1; 
%       insig = x(Ni:Nf);
%       p = 180;
%       Get_ESPRIT_analysis(insig,p,N,fs);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20151204_Antoine_lecture
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 05/12/2015
% Last update on: 05/12/2015 
% Last use on   : 06/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    p = 180;
end

X0 = il_get_X0(insig,p,'X0'); % Laroche1993, Eq. 2. Matrix (N-p) x p
X1 = il_get_X0(insig,p,'X1'); % Laroche1993, Eq. 3. Matrix (N-p) x p
X0plus = pinv(X0); % Laroche1993, Step (2)

% Laroche1993, Step (3):
Matr = X0plus*X1;  
lambda = eig(Matr); 
    
% Laroche1993, Step (4):
idxposim = find(imag(lambda)>0); % finds lambda with imaginary part > 0
    
L = length(idxposim);

%%% Frequency and damping factors:
fi = ( atan2(imag(lambda(idxposim)),real(lambda(idxposim))) )/(2*pi);
alpha_i = abs(-log(abs(lambda(idxposim))));

%%% Form matrix T
T = nan(M,2*L);

for n = 1:M
    nu_i(n,:)   = transpose(exp(-n*alpha_i) .* cos(2*pi*n*fi));
    beta_i(n,:) = transpose(exp(-n*alpha_i) .* sin(2*pi*n*fi));
end

idxs = 1:2:2*L;
T(1,idxs) = 1;

idxs = 2:2:2*L;
T(1,idxs) = 0;

for n = 1:L
    T(2:end,2*n-1) = nu_i(2:end,n);
    T(2:end,2*n  ) = beta_i(2:end,n);
end

Ttra = T';

MultM = pinv(Ttra*T) * Ttra*insig(1:M);
ai = MultM(1:2:end);
bi = MultM(2:2:end);

Ai = sqrt(ai.^2+bi.^2);
Phi_i = tan(bi./ai);

%%% Obtained parameters: Ai, alpha_i, fi, Phi_i

n = transpose(1:fs/2);
outsig = zeros(size(n));

[xx maxidx] = max(Ai);

[fi fi_idx] = sort(fi);
Ai      = Ai(fi_idx);
alpha_i = alpha_i(fi_idx);
Phi_i   = Phi_i(fi_idx);

Fi      = fi*fs;
sigma_i = alpha_i*fs;

if nargout == 0
    figure;
    stem(Fi,Ai); grid on
    xlabel('Frequency [Hz]')
end

for i = 1:L %maxidx %L

    fac1 = Ai(i);
    fac2 = exp(-alpha_i(i)*n);
    fac3 = cos(2*pi*n*fi(i)+Phi_i(i));
    outsig = outsig+fac1*fac2.*fac3;
end

if nargout == 0
    figure;
    subplot(2,1,1)
    plot(outsig,'r');
    title('Synthesised signal')
    
    subplot(2,1,2)
    x2plot = insig;
    plot(x2plot);
    title('Input signal')

    sound(outsig,fs)
    pause(1.2*length(outsig)/fs)
    sound(x2plot,fs)
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
