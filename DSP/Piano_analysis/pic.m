function XOUT  = pic(XIN,SEUIL,MIN)   

%usage
%    function XOUT  = pic(XIN,SEUIL,MIN)
%RETOURNE 
% XPEAK     : LE MODULE DU SPECTRE DE FOURIER APRES DETECTION DE PICS
%             LE MODULE EST EN LINEAIRE ....
%ARGUMENT
% XIN        MATRICE DES FFT EFFECTUEES SUR CHAQUE FENETRE
% SEUIL      SEUIL DE RECHERCHE (<0)
% MIN        VALEUR MINIMALE DE L"AMPLITUDE DES 
%            PICS SELECTIONNABLES, MIN EST LINEAIRE
%            ET SA VALEUR EST CALCULEE GRACE AUX 
%            SPECTRES PRECEDENTS
%            MIN=max du module du spectre precedent - X dB 
%            X dB etant > SEUIL generalement 
% remarque un seul critere retenu

N = length(XIN);   

if nargin == 2
    MIN = 0;
end

% DECALAGE VERS LA DROITE  
XMOINS1 = [0,XIN(1:N-1)]  ;
% DECALAGE VERS LA GAUCHE  
XPLUS1  = [XIN(2:N),0]  ;
CRITERE1_2 = XIN>=XPLUS1 & XIN>XMOINS1  ;
XOUT = XIN .* CRITERE1_2  ;

SEUIL  = ones(1,N) .* (10^( ( 20*log10(max(XOUT)) + SEUIL ) /20) ) ; 
CRITERE_SEUIL = XOUT>SEUIL   ; 
XOUT = XOUT .* CRITERE_SEUIL    ;  


% SEUIL DE DISCRIMINATION GLOBAL
CRITERE_GLOB = XOUT>MIN   ; 
XOUT = XOUT .* CRITERE_GLOB    ;  

% SUPPRESSION DE LA RAIE A L'INDICE 1 ET A L'INDICE N 
XOUT(1) = 0     ;
XOUT(N) = 0      ;

