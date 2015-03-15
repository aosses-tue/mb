% Script RoughBody
% 
% Created by: Dik Hermes
% Received in: August 2014
% Additional comments by: Alejandro Osses
% Last edit on: -
% Last used on: 15/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window	= blackman(N);
AmpCal	= db2amp(80)*2/(N*mean(window));	% Calibration between wav-level and loudness-level (assuming blackman window and FFT will follow)
window	= AmpCal*(window)';		% -20 dBFS	<-> 60 dB SPL
Chno	= 47;
Cal		= 0.25;
N2		= N/2;
q		= 1:1:N;
qb		= N0:1:Ntop;
freqs	= (qb+1)*Fs/N;
hBPi	= zeros(Chno,N);
hBPrms	= zeros(1,Chno);
mdept	= zeros(1,Chno);
ki		= zeros(1,Chno-2);
ri		= zeros(1,Chno);

% Calculate Excitation Patterns
TempIn	= window.*(InputSig');
TempIn	= a0.*fft(TempIn);
Lg		= abs(TempIn(qb));
LdB		= amp2db(Lg);
whichL	= find(LdB>MinExcdB);
sizL	= length(whichL);

% steepness of slopes (Terhardt)
S1		= -27;			
S2		= zeros(1,sizL);
for w	= 1:1:sizL;
   steep= -24-(230/freqs((w)))+(0.2*LdB(whichL(w)));	%Steepness of upper slope [dB/Bark] in accordance with Terhardt
   if steep<0
      S2(w)	= steep;
   end
end

whichZ      = zeros(2,sizL);
qd          = 1:1:sizL;
whichZ(1,:)	= floor(2*Barkno(whichL(qd)+N01));
whichZ(2,:)	= ceil(2*Barkno(whichL(qd)+N01));

ExcAmp = zeros(sizL,47);
Slopes = zeros(sizL,47);

for k=1:1:sizL
    
    Ltmp = LdB(whichL(k));
	Btmp = Barkno(whichL(k)+N01);
    for l = 1:1:whichZ(1,k)
        Stemp =	(S1*(Btmp-(l*0.5)))+Ltmp;
        if Stemp>MinBf(l)
            Slopes(k,l)=db2amp(Stemp);
        end
    end
    for l=whichZ(2,k):1:47
        Stemp	=	(S2(k)*((l*0.5)-Btmp))+Ltmp;
        if Stemp>MinBf(l)
            Slopes(k,l)=db2amp(Stemp);
        end
    end
    
end

for k=1:1:47
    etmp=zeros(1,N);
    for l=1:1:sizL
        N1tmp=whichL(l);
        N2tmp=N1tmp+N01;
        if whichZ(1,l)==k
            ExcAmp(N1tmp,k)	=	1;   
        elseif whichZ(2,l)==k
            ExcAmp(N1tmp,k)	=	1;   
        elseif whichZ(2,l)>k 
            ExcAmp(N1tmp,k)	=	Slopes(l,k+1)/Lg(N1tmp);
        else 
            ExcAmp(N1tmp,k) =	Slopes(l,k-1)/Lg(N1tmp);      
        end
            etmp(N2tmp)		=	ExcAmp(N1tmp,k)*TempIn(N2tmp);      
    end    
    ei(k,:)		=	N*real(ifft(etmp));
	etmp		=	abs(ei(k,:));
    h0(k)		=	mean(etmp);
    Fei(k,:)	=	fft(etmp-h0(k));
    hBPi(k,:)	=	2*real(ifft(Fei(k,:).*Hweight(k,:)));	
	hBPrms(k)	=	rms(hBPi(k,:));   
    if h0(k)>0
        mdept(k)=	hBPrms(k)/h0(k);
        if mdept(k)>1
            mdept(k)=1;
        end
    else
        mdept(k)=0;
    end

    Roughness('Prog',47*ProgMax,(ProgCur-1)*47+k); %show progress bar

end

% find cross-correlation coefficients
for k=1:1:45
   cfac	=	cov(hBPi(k,:),hBPi(k+2,:));
   den	=	diag(cfac);   
   den	=	sqrt(den*den');
   if den(2,1)>0
		ki(k)	=	cfac(2,1)/den(2,1);
   else
      ki(k)	=	0;
   end   
end

% Calculate specific roughness ri and total roughness R
ri(1)	=	(gzi(1)*mdept(1)*ki(1))^2;
ri(2)	=	(gzi(2)*mdept(2)*ki(2))^2;
for k=3:1:45
   ri(k)=	(gzi(k)*mdept(k)*ki(k-2)*ki(k))^2;
end
ri(46)	=	(gzi(46)*mdept(46)*ki(44))^2;
ri(47)	=	(gzi(47)*mdept(47)*ki(45))^2;

R			=	Cal*sum(ri);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%