function [A, F]= freq_xpeak(X, Xpeak, fe)   

% Usage
%
% [A, F] = freq_xpeak (X, Xpeak, fe)
%
% Retourne :
%
%   A         Amplitude corrigee des raies
%   F         Frequence corrigee des raies
%
% Les arguments :
%
%   X         Spectre initial  
%   Xpeak     Spectre apres detection de pics 
%   fe        Frequence d'echantillonnage  

% Passage au log

Xlog = 20*log10(X);
 
for i=1:length(X), 
% Attention! compensation de 1 au niveau de l'indice pour
% retrouver la frequence reelle    
  
  if Xpeak(i) > 0
      x = [(i-2) (i-1) i]';
      y = [Xlog(i-1) Xlog(i) Xlog(i+1)]';   
      c = polyfit(x, y, 2);   
      F(i) = -c(2)/(2*c(1));
      A(i) = 10^(polyval(c, F(i))/20);
  else
      F(i) = 0;
      A(i) = 0; 
  end
end 

% Nombre de points pour le calcul de la FFT

nb_points = 2*(length(X) - 1) ;
 
% Resolution frequentielle
  
r = fe/nb_points;

% Frequence en Hertz
 
F = F*r; 

