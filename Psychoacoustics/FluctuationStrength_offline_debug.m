function dataOut = FluctuationStrength_offline_debug(dataIn, Fs, N, bDebug)
% function dataOut = FluctuationStrength_offline_debug(dataIn, Fs, N, bDebug)
%
% 1. Description:
%       Off-line implementation of the roughness algorithm.
%       Changes:
%           amp2db replaced by To_dB
%           db2amp replaced by From_dB
%           private rms renamed to dw_rms
%
% author : Matt Flax <flatmax @ http://www.flatmax.org> : Matt Flax is flatmax
% March 2006 : For the psySoundPro project
%
% revised : Farhan Rizwi
%           July '07
%           Reformatted, copied and vectorised code from InitAll
%           and Hweights into the function space below.  This
%           allows us to use nested functions effeciently.
%
% contact for the original source code :
% http://home.tm.tue.nl/dhermes/
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 10/11/2014
% Last update on: 10/11/2014 % Update this date manually
% Last use on   : 10/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    bDebug = 0;
end

if nargin < 3
    N = 8192;
end

if ~(Fs == 44100 | Fs == 40960 | Fs == 48000)
  warning(['Incorrect sample rate for this roughness algorithm. Please ' ...
           're-sample original file to be Fs=44100,40960 or 48000 Hz']);
end

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%

Bark = Get_psyparams('Bark');
    
Bark2   = [sort([Bark(:,2);Bark(:,3)]),sort([Bark(:,1);Bark(:,4)])];

N0      = round(20*N/Fs)+1;
N01     = N0-1;
N50     = round(50*N/Fs)-N0+1;
N2      = N/2+1;
Ntop	= min( round(20000*N/Fs)+1, round(N/2) ); % mod by AO
Ntop2	= Ntop-N0+1;
dFs	= Fs/N;

% Make list with Barknumber of each frequency bin
Barkno      = zeros(1,N2);
f           = N0:1:Ntop;
Barkno(f)   = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

% Make list of frequency bins closest to Cf's
Cf = ones(2,24);
for a=1:1:24
    
    Cftmp = Bark(a+1,2);
    
    if Cftmp < Fs/2
        Cf(1,a)=round(Bark((a+1),2)*N/Fs)+1-N0;
        Cf(2,a) = Cftmp;
    end
    
end
%Make list of frequency bins closest to Critical Band Border frequencies
Bf = ones(2,24);
Bf(1,1)=round(Bark(1,3)*N/Fs);
for a=1:1:24
    
    Ctmp = Bark((a+1),3);
    
    if Ctmp < Fs/2
        Bf(1,a+1)=round(Ctmp*N/Fs)+1-N0;
        Bf(2,a)=Bf(1,a)-1;
    end
end

idx = min( max( find(Bark(:,3) <= 20000)) ,max( find(Bark(:,3) < Fs/2)));
Bf(2,25)=round(Bark(idx,3)*N/Fs)+1-N0;

%Make list of minimum excitation (Hearing Treshold)
HTres = Get_psyparams('HTres');

k = (N0:1:Ntop);
MinExcdB = interp1(HTres(:,1),HTres(:,2),Barkno(k));
  
% Initialize constants and variables
zb    = sort([Bf(1,:),Cf(1,:)]);

idxtmp = find(zb==0 | zb==1);
zb(idxtmp) = [];

MinBf = MinExcdB(zb);

Chno    = min(idx*2-1,47);
ei      = zeros(Chno,N);
Fei     = zeros(Chno,N);

h0     = zeros(1,Chno);
k      = 1:1:Chno;
    
% calculate a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 9; % similar setting than Roughness, see Sontacchi1999, pp 76 of 115
switch method
    case 1
        kg = 1;     s = 1;      p = 1;      qg = 0;
    case 2
        kg = 1;     s = 1;      p = 1.5;    qg = 0;
    case 3
        kg = 1;     s = 1;      p = 2;      qg = 0;
    case 4
        kg = 1;     s = 1;      p = 1;      qg = 0.07;
    case 5
        kg = 1;     s = 1;      p = 1.5;    qg = 0.07;
    case 6
        kg = 1;     s = 1;      p = 2;      qg = 0.07;
    case 7
        kg = 2;     s = 1;      p = 1;      qg = 0;
    case 8
        kg = 2;     s = 1;      p = 1.5;    qg = 0;
    case 9
        kg = 2;     s = 1;      p = 2;      qg = 0;
    case 10
        kg = 2;     s = 1;      p = 1;      qg = 0.07;
    case 11
        kg = 2;     s = 1;      p = 1.5;    qg = 0.07;
    case 12
        kg = 2;     s = 1;      p = 2;      qg = 0.07;
end

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

% BEGIN Hweights:

Htmp = Get_Hweight_roughness(N,Fs);
Hweight = repmat( Htmp(17,:),Chno,1);
% Hweight = Get_Hweight_roughness(N,Fs);

% END Hweights 
%%%%%%%%%%%%%%%%

% BEGIN: FluctBody

bBlackman = 1;

Window = blackman(N, 'periodic') .* 1.8119;
dBcorr = 80+10.72; % originally = 80

if bBlackman
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 
else
    warning('Blackman by-passed')
    AmpCal = From_dB(dBcorr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration between wav-level and loudness-level (assuming
% blackman window and FFT will follow)
    
% Chno	=	47;
Cal	 	=	0.25;
N2		=	N/2;
q		=	1:1:N;
qb		=	N0:1:Ntop;
freqs	=	(qb+1)*Fs/N;
hBPi	=	zeros(Chno,N);
hBPrms	=	zeros(1,Chno);
mdept	=	zeros(1,Chno);
ki		=	zeros(1,Chno-2);
fi		=	zeros(1,Chno);

% Calculate Excitation Patterns
TempIn  = dataIn*AmpCal;
[rt,ct] = size(TempIn);
[r,c]   = size(a0);
if rt  ~= r; TempIn=TempIn'; end
    
TempIn	= a0.*fft(TempIn);
Lg		= abs(TempIn(qb));
LdB		= To_dB(Lg);
whichL	= find(LdB>MinExcdB);
sizL	= length(whichL);

% steepness of slopes (Terhardt)
S1 = -27;
S2 = zeros(1,sizL);

for w = 1:1:sizL;
    % Steepness of upper slope [dB/Bark] in accordance with Terhardt
    steep = -24-(230/freqs(w))+(0.2*LdB(whichL(w)));
    if steep < 0
        S2(w) = steep;
    end
end

whichZ	= zeros(2,sizL);
qd		= 1:1:sizL;
whichZ(1,:)	= floor(2*Barkno(whichL(qd)+N01));
whichZ(2,:)	= ceil(2*Barkno(whichL(qd)+N01));

ExcAmp = zeros(sizL,Chno);
Slopes = zeros(sizL,Chno);

for k=1:1:sizL
    Ltmp = LdB(whichL(k));
    Btmp = Barkno(whichL(k)+N01);

    for l = 1:1:whichZ(1,k)
        Stemp = (S1*(Btmp-(l*0.5)))+Ltmp;
        if Stemp>MinBf(l)
          Slopes(k,l)=From_dB(Stemp);
        end
      end
      for l = whichZ(2,k):1:Chno
        Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp;
        if Stemp>MinBf(l)
          Slopes(k,l)=From_dB(Stemp);
        end
      end
end

for k=1:1:Chno
    
    etmp = zeros(1,N);
    
    for l=1:1:sizL % each level above threshold
        
        N1tmp = whichL(l);
        N2tmp = N1tmp + N01;
        if (whichZ(1,l) == k)
            ExcAmp(N1tmp, k) = 1;
        elseif (whichZ(2,l) == k)
            ExcAmp(N1tmp, k) = 1;
        elseif (whichZ(2,l) > k)
            ExcAmp(N1tmp,k) = Slopes(l,k+1)/Lg(N1tmp);
        else
            ExcAmp(N1tmp,k) = Slopes(l,k-1)/Lg(N1tmp);
        end
        etmp(N2tmp) = ExcAmp(N1tmp,k)*TempIn(N2tmp); 
    end
    
    if k == 17
        if bDebug
            figure(2); 
            plot(freqs,20*log10(abs(etmp(qb))))
            hold on; grid on
            xlabel('Frequency [Hz]')
            ylabel('Level [dB]')
            title('Spectral components: Terhardt''s model')
            xlim([900 1100])
        end
    end
    
    % etmp_fd  - excitation pattern in frequency domain
    % etmp  - excitation patterns in time domain (after L242)
    % ei    - excitation patterns in time domain
    % Fei   - envelope in frequency domain, as function of fmod
    % h0    - DC component, average of half-wave rectified signal
    
    etmp_fd(k,:)   = etmp;
    ei(k,:)     = N*real(ifft(etmp));
    etmp_td(k,:)= abs(ei(k,:));
    h0(k)       = mean(etmp_td(k,:));
    Fei(k,:)	= fft( etmp_td(k,:)-h0(k) ); % changes the phase but not the amplitude
    
    if bDebug
        Fsnew       = Fs/2;
        Nnew        = N/2;
        eir         = resample(ei(k,:),Fsnew,Fs); % eir at Fs = 22050
        etmpr       = abs(eir);
        h0r(k)      = mean(etmpr);
        
        % Feir(k,:)	= fft( etmpr-h0r(k) ,Nnew); % TEST
        Feir(k,:)	= fft( etmpr-h0r(k) ,N); 
        Hweightr(k,:) = Hweight(k,:); % resample( Hweight(k,:),Fsnew,Fs );
    
        % hBPi  - band-pass filtered envelopes in time domain
    
        [b,a] = butter(8,1000/Fs,'low');
        [bred,ared] = butter(8,(Fsnew/2)/Fsnew,'low');
        eirFS = filter(bred,ared,eir);
        h0rFS(k) = mean(abs(eirFS));
    end
    
    hBPi(k,:)	= 2*real(  ifft( Fei(k,:).*Hweight(k,:) )  );
    hBPrms(k)	= dw_rms(hBPi(k,:));
    
    if bDebug
        hBPi2 = filter(b,a,hBPi(k,:));
        hBPrms2(k)	= dw_rms(hBPi2);
        
        hBPi3(k,:)	= 2*real(  ifft( Feir(k,:).*Hweightr(k,:) ,N)  );
        hBPi3(k,:)  = filter(bred,ared,hBPi3(k,:));
        hBPrms3(k)	= dw_rms(hBPi3(k,:));
    end
    
    if k == 17
        if bDebug
            
            fidx_min = min(qb);
            figure(3); 
            
            Fei_dB = 20*log10(abs(Fei(k,:)));
            plot(freqs, Fei_dB(qb),'LineWidth',1), hold on
            
            opts.fs = Fs;
            [ytmp tmpYdB,ftmp] = freqfft(hBPi(k,:)',N/2,opts);
            
            plot(ftmp(fidx_min:end),tmpYdB(fidx_min:end),'r')
            grid on
            xlabel('Frequency [Hz]')
            ylabel('Level [dB]')
            
            opts.fs = Fs;
            [ytmp tmpYdB2,ftmp] = freqfft(hBPi2',N/2,opts);
            plot(ftmp(fidx_min:end),tmpYdB2(fidx_min:end),'m--')
            
            opts.fs = Fsnew;
            [ytmp tmpYdB3,ftmp3] = freqfft(hBPi3(k,:)',N/2,opts);
            plot(ftmp3,tmpYdB3,'k--','LineWidth',1)
            
            legend('Fei','hBPi','hBPi with LPF','hBPi with Fs_r_e_d + LPF')
        end
    end
    
    if h0(k)>0
        mdept(k) = hBPrms(k)/h0(k);
        if bDebug
            mdept3(k) = hBPrms3(k)/h0r(k);
        end
        if mdept(k)>1
            mdept(k)=1;
        end
        if bDebug
            mdept3(k) = 1;
        end
    else
        mdept(k)=0;
    end
    
end

if bDebug
    
    % params for figures 4 and 5
    idx2plot = 14:22; 
        
    % params just for figure 4
    idxq = find(freqs<1200 & freqs>800);
    L = length(idxq);
    freqsm = repmat(freqs(idxq),length(idx2plot),1);
    BandNumber = repmat( ((idx2plot)'), 1, L);
    
    figure(4);
    mesh( freqsm,BandNumber,20*log10(abs(etmp_fd(idx2plot,idxq))) )
    xlabel('Frequency [Hz]')
    ylabel('Band number')
    zlabel('Level [dB]')
    title('Excitation patterns')
    
    % params just for figure 5
    BandNumber = repmat( ((idx2plot)'), 1, N);
    t = repmat( (1:N)/Fs ,length(idx2plot),1);
    
    figure(5); 
    mesh( t,BandNumber,hBPi(idx2plot,:) )
    xlabel('Time [s]')
    ylabel('Band number')
    title('Band-pass signals')
    zlim([-4000 4000])
    
end

% find cross-correlation coefficients
    for k=1:1:Chno-2
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
    fi(1)	=	(h0(1)^qg)*(mdept(1)^p)*(ki(1)^kg);
    fi(2)	=	(h0(2)^qg)*(mdept(2)^p)*(ki(2)^kg);
    for k = 3:1:Chno-2
      fi(k)	=	(h0(k)^qg)*(mdept(k)^p)*( ki(k-2)*ki(k) )^kg;
    end
    fi(Chno-1)	=	(h0(Chno-1)^qg)*(mdept(Chno-1)^p)*(ki(Chno-3)^kg);
    fi(Chno  )	=	(h0(Chno  )^qg)*(mdept(Chno  )^p)*(ki(Chno-2)^kg);
    FS		=	Cal*( sum(fi) )^(1/s);
    
    SPL = mean(rms(dataIn));
    if SPL > 0
      SPL = To_dB(SPL)+dBcorr+3; % -20 dBFS <--> 60 dB SPL
    else
      SPL = -400;
    end
    
    % Create a cell array to return
    dataOut{1} = FS;
    dataOut{2} = fi;
    dataOut{3} = SPL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: RoughBody

end
