function [outsig fi ai alpha_i Phi_i] = Get_ESPRIT_analysis(insig,p,N,fs)
% function [outsig Fi Ai Alpha_i Phi_i] = Get_ESPRIT_analysis(insig,p,N,fs)
%
% 1. Description:
%       outsig ois a synthesised signal by approximating the insig signal
%       using the pencil matrix method.
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/12/2015
% Last update on: 05/12/2015 
% Last use on   : 05/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    p = 180;
end

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
alphai = abs(-log10(abs(lambda(idxposim))));
Alphai = alphai*fs;
    
M = N; % this should not be necessarily

%%% Form matrix T
T = nan(M,2*L);

for n = 1:M
    nu_i(n,:)   = transpose(exp(-n*alphai) .* cos(2*pi*n*fi));
    beta_i(n,:) = transpose(exp(-n*alphai) .* sin(2*pi*n*fi));
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

if nargout == 0
    figure;
    stem(Fi,Ai); grid on
end

n = transpose(1:fs/2);
outsig = zeros(size(n));

[xx maxidx] = max(Ai)
for i = 1:L %maxidx %L

    fac1 = Ai(i);
    fac2 = exp(-alphai(i)*n);
    fac3 = cos(2*pi*n*fi(i)+Phi_i(i));
    outsig = outsig+fac1*fac2.*fac3;
end

if nargout == 0
    figure;
    subplot(2,1,1)
    plot(outsig,'r');

    subplot(2,1,2)
    x2plot = insig;
    plot(x2plot);

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
