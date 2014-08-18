function y = ch_testnoise(l,fu,fo,t,option)
% function y = ch_testnoise(l,fu,fo,t,option)
%
% y = ch_testnoise(l,fu,fo,t,option)
% berechnet Weisses Rauschen mit Pegel l (in dB), unterer Grenzfrequenz fu (in Hz), 
% oberer Grenzfrequenz fo (in Hz)und Dauer t (in s).
% Achtung! muss geaendert werden (s. medavnoise) asymmetrische Filter
% falls 'g' als weiterer Parameter übergeben wird, wird das Rauschen gaußmoduliert (20 ms)
 
Fs= 44100;
Fn= Fs/2;
noise=randn(t*Fs,1);          % gaussverteiltes Rauschen mean=0, std=1 und Dauer t

% Bandbegrenzung
fm = sqrt(fu*fo);
Wpl= fo/Fn; % passband normalised freq L
Wph= fu/Fn; % passband normalised freq H
Rs = 60;
Rp = 1;

zm=tonheit(fm);
zo=60/27+zm;
zu=zm-60/27;
Wsl=invtonheit(zo)/Fn; % stopband freq
Wsh=invtonheit(zu)/Fn; % stopband freq

%Tiefpassfilterung
[n,Wn]=cheb1ord(Wpl,Wsl,Rp,Rs);
[b,a]=cheby1(n,Rp,Wn);
noise=filter(b,a,noise);

%Hochpassfilterung
[n,Wn]=cheb1ord(Wph,Wsh,Rp,Rs);
[b,a]=cheby1(n,Rp,Wn,'high');
noise=filter(b,a,noise);

rms=true_rms(noise);
x=1;
while (rms<l-0.5) | (rms>l+0.5)   % gewuenschter Pegel wird eingestellt
   if rms<l-0.5
      x=x+x/2;
   else
      x=x-x/2;
   end
   sig=x.*noise;
   rms=true_rms(sig);
end   
y=x*noise;

if (nargin == 5)         % Gaussmodulation (falls gewuenscht)
   if strcmp(option,'g')        
      y=gaussmod(y,20);
   else
      error('falscher Parameter!')
   end
end
