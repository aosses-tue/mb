function [O f] = spect(insig,First,siz,fs,cc,normal,N)
% function [O f] = spect(insig,First,siz,fs,cc,normal,N)
%
% 1. Description:
%     spect(z,First,siz,Fe) trace la FFT du signal z
%     en commencant a l'echantillon begin, et en utilisant siz echantillons
%     Une fenetre de Hann est pre-appliquee, et le
%     spectre de puissance est retourne, en log avec un plancher de -100dB.
%     spec(z,begin,siz,fs,'--') trace en pointilles au lieu des traits pleins.
%     spec(z,begin,siz,fs,'-',0) effectue la meme operation sans normalisation
%     pratique pour comparer deux signaux du point de vue puissance.
%     enfin spec(z,begin,siz,Fe,'-',1,N) fait une fft en utilisant siz echantillons
%     mais avec un zero padding jusqu'a N.
%     O = spect(z,begin,siz,Fe) retourne le vecteur spectre.
%MAX
% 2. Stand-alone example:
%     file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%     [insig fs] = Wavread(file);
%     N = 2^11;
%     Ni  = 3191; % max of the waveform (manually computed)
%     Nf  = Ni+N-1; 
%     figure;
%     spect(insig,Ni,N,fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin == 4) 
    normal = 1 ; 
    cc = '-b'; 
    N = siz;
end 

if(nargin == 5) 
    normal = 1; 
    N = siz;
end 

if(nargin == 6) 
    N = siz;
end 

if(N < siz) 
    error('Attention N est trop petit') ;
end

x = insig(First:First+siz-1);
%h = blackman(siz)';
h = hanning(siz);
x = x .* h ;
X = abs(fft(x,N));
X = 20 * (log10(X) - normal * log10(max(X)));
X = max(X,-100);
if(normal == 1)  
    axis([0 fs*.5 -100 0]);
end

f = fs*(0:1/N:.5-1/N);
plot(f,X(1:N/2),cc);

if(normal == 1) 
    axis;
end 

xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');

O = X(1:N/2);
