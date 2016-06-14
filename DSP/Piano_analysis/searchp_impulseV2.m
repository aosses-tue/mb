function [AdB, F] = searchp_impulseV2(insig, fs, First, End, SEUIL)
% function [AdB, F] = searchp_impulseV2(insig, fs, First, End, SEUIL)
% 
% 1. Description:
%       Fonction de recherche de pics frequentiels contenus dans un signal
%       Elabore pour l'analyse des impulsions - projet PAPA
% Auteur :  Antoine Chaigne - mai 2015 -
% usage :
%   [AdB, F] = searchp_impulseV2(sig, fe, First, End, SEUIL)
%
% Les arguments :
%
% insig   signal a analyser 
% fs      Frequence d'echantillonnage (en Hz)
% First   Premier echantillon (par defaut : First = 1)
% End     Dernier echantillon (par defaut : End = 2048)
% SEUIL   Valeur (en dB) pour la detection de pics (typ. = -40)
% AdB     Amplitude (en dB) des pics detectes
% F       Frequence (en Hz) des pics detectes
%
% 2. Stand-alone example:
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [insig fs] = Wavread(file);
%       N = 2^11;
%       Ni  = 3191; % max of the waveform (manually computed)
%       Nf  = Ni+N-1; 
%       figure;
%       searchp_impulseV2(insig,fs,Ni,Nf);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    SEUIL = -40; % -40 dB
end

if (nargin == 3)
    End = First + 2^11;
elseif (nargin == 2)
    First = 1;
    End = First + 2^11;
end

if size(insig,2) == 1
    insig = transpose(insig);
    warning('A row vector is expected as input signal');
end
%Calcul de la valeur efficace du signal
sigeff=sqrt(sum(insig.*insig))/sqrt(length(insig));
%Normalisation
insig=insig/sigeff;
%FFT avec fenetre de Hanning 
XIN = module_han(insig(First:End));
%Detection de pics
XOUT  = pic(XIN,SEUIL);
%Transformation des pts FFT en frequence
[A,F] = freq_xpeak(XIN, XOUT, fs);
% Elimination des valeurs nulles
A=A(A~=0);
F=F(F~=0);

%Calcul de l'amplitude des pics detectes en dB
AdB = 20*log10(A/max(A));
if nargout == 0
    %Representation graphique A(F)
    % plot_xpeak_ONL(A, F, fs, SEUIL);
    plot_xpeak(A, F, fs, SEUIL);
    spectre=[F' A'];
end

if nargout == 0
    plot(F,AdB,'*');
end
