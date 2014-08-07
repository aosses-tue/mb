function y = Create_test_noise(l,fu,fo,t)
% function y = Create_test_noise(l,fu,fo,t)
%
% 1. Description:
%       Inspired in Chalupper's function 'testnoise.m'
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       sig = Create_test_noise(70,800,1200,1);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs  = 44100;
fn  = fs/2;             % Nyquist freq
noise=randn(t*fs,1);    % gaußverteiltes Rauschen mean=0, std=1 und Dauer t

% Bandbegrenzung
fm = sqrt(fu*fo);
Wpl= fo/fn; % passband normalised freq L
Wph= fu/fn; % passband normalised freq H
Rs = 60;
Rp = 1;

try
    zm =tonheit(fm); % Original Chalupper's function
    zu =zm-60/27;
    zo =60/27+zm;
catch
    zm = hz2mel(fm); 
    zu =zm-hz2mel(fu);
    zo =hz2mel(fo)+zm;
end
    
try
    Wsl=invtonheit(zo)/fn; % stopband freq
    Wsh=invtonheit(zu)/fn; % stopband freq
catch
    Wsl=mel2hz(zo)/fn; % stopband freq
    Wsh=mel2hz(zu)/fn; % stopband freq
end
 
if Wsl > 1
    Wsl = 1-1e-4;
end

if Wsl < 0
    Wsl = 0+1e-4;
end

%Tiefpaßfilterung
[n,Wn]=cheb1ord(Wpl,Wsl,Rp,Rs);
[b,a]=cheby1(n,Rp,Wn);
noise=filter(b,a,noise);

%Hochpaßfilterung
[n,Wn]=cheb1ord(Wph,Wsh,Rp,Rs);
[b,a]=cheby1(n,Rp,Wn,'high');
noise=filter(b,a,noise);

try
    rms=true_rms(noise);
catch
    rms = dbspl(noise);
end

x=1;
while (rms<l-0.5) | (rms>l+0.5)   % gewünschter Pegel wird eingestellt
    % AO: maybe quicker if current 'calibration' function is used instead
   if rms<l-0.5
      x=x+x/2;
   else
      x=x-x/2;
   end
   sig=x.*noise;
   try
        rms=true_rms(sig);
   catch
       rms = dbspl(sig);
   end
end   
y=x*noise;

if (nargin == 5)         % Gaußmodulation (falls gewünscht)
   if strcmp(option,'g')        
      y=gaussmod(y,20);
   else
      error('falscher Parameter!')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
