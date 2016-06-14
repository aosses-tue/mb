function [Yf] = stft(strwav,N,Nf,ch,fmax)
% function [Yf] = stft(strwav,N,Nf,ch,fmax)
%
% [Yf] = stft(strwav,N,Nf,ch,fmax)
% strwav : nom du fichier .wav à analyser
% N : nombre de fenêtres sur lesquelles on veut calculer la FFT
% Nf : taille de la fenêtre en nombre d'échantillons
% ch : numero de la piste à analyser dans le fichier .wav
% fmax (optionnel) : fréquence maximale de l'affichage : on obtient un graphe sur
% l'intervalle [0 fmax], par défaut fmax=2000
%
% Le recouvrement des différentes fenêtres se gère en ajustant leur nombre
% et leur taille par rapport à la taille totale du signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('Still to analyse by AO')
Fe=51200; %Frequence d'echantillonnage

if (nargin == 4)
    fmax=2000;
end

X=audioread(strwav);
X=X(:,ch)';

s=size(X,2);
p=round((s-Nf)/(N-1));
Yf=0;
for k=0:N-2
    Y=[zeros(1,k*p), hann(Nf)'.*X((k*p+1):(k*p+Nf)), zeros(1,s-(k*p+Nf))];
    Yf=Yf+abs(fft(Y));
end

Y=[zeros(1,s-Nf),hann(Nf)'.*X((s-Nf+1):s)];
Yf=Yf+abs(fft(Y));
M=max(Yf);
plot(Fe/s*[1:round(fmax*s/Fe)],Yf(1:round(fmax*s/Fe))/M);

end
