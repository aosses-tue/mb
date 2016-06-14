function X = module_han(x)

%usage
%    function X = module_han(x)
% RETOURNE LE MODULE DE LA FFT SUR N/2 + 1 POINTS  
% FENETRE DE PONDERATION => HANNING   

x = x .* hanning(length(x))'  ;
f=fft(x)  ;
module = abs(f);

X = module(1:(length(f)/2)+1)  ;

