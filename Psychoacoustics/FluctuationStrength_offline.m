function dataOut = FluctuationStrength_offline(dataIn, Fs, N, bDebug)
% function dataOut = FluctuationStrength_offline(dataIn, Fs, N, bDebug)
%
% 1. Description:
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
  error(['Incorrect sample rate for this roughness algorithm. Please ' ...
         're-sample original file to be Fs=44100,40960 or 48000 ' ...
         'Hz']);
end

Fs_red  = 20; % in Hz, Added by AO

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%
Bark = Get_psyparams('Bark');

Bark2	= [sort([Bark(:,2);Bark(:,3)]),sort([Bark(:,1);Bark(:,4)])];
fi      = 1; % fi = 20 for roughness
N0      = round(fi*N/Fs)+1;
N01     = N0-1;
N50     = round(50*N/Fs)-N0+1;
N2      = N/2+1;
Ntop	= round(20000*N/Fs)+1; % max f = 20 kHz
% Ntop	= N;
Ntop2	= Ntop-N0+1;
dFs     = Fs/N;
dFs_red = Fs_red/N;

% Make list with Barknumber of each frequency bin
Barkno      = zeros(1,N2);
f           = N0:1:Ntop;
Barkno(f)   = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

f_red = (1:N)/N * Fs_red;

% Make list of frequency bins closest to Cf's
Cf = ones(2,24);
for a=1:1:24
  Cf(1,a)=round(Bark((a+1),2)*N/Fs)+1-N0;
  Cf(2,a)=Bark(a+1,2);
end
%Make list of frequency bins closest to Critical Band Border frequencies
Bf = ones(2,24);
Bf(1,1)=round(Bark(1,3)*N/Fs);
for a=1:1:24
  Bf(1,a+1)=round(Bark((a+1),3)*N/Fs)+1-N0;
  Bf(2,a)=Bf(1,a)-1;
end
Bf(2,25)=round(Bark((25),3)*N/Fs)+1-N0;

%Make list of minimum excitation (Hearing Treshold)
HTres= Get_psyparams('HTres');

k = (N0:1:Ntop);
MinExcdB = interp1(HTres(:,1),HTres(:,2),Barkno(k));
  
% Initialize constants and variables
zi    = 0.5:0.5:23.5;
zb    = sort([Bf(1,:),Cf(1,:)]);
MinBf = MinExcdB(zb);
ei    = zeros(47,N);
Fei   = zeros(47,N);

h0     = zeros(1,47);
k      = 1:1:47;
    
% calculate a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

method = 1; % arbitrarily, see Sontacchi1999, pp 76 of 115
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN Hweights, by AO
% weights for freq. bins < N/2
Hweight	= zeros(1,N);

DCbins	= 0;

H2 = [	0       0
        0.2441  0
        0.7324	0.45
        1.7090	0.82
        1.9531  0.9
        3.1738	1
        12.9395 1
        14.4043 0.91
        15.3809	0.55
        15.8691 0.25  
        16      0];     % obtained from Sontacchi1999, Abbildung 5.4, understanding
                        % that 'BINS' are converted to frequency in such a way that
                        % bin 16th is approximately 16 Hz
                        
fH2 = H2(:,1); 

last	= floor((16/Fs_red)*N) ;
k       = DCbins+1:1:last;
f       = (k-1)*Fs_red/N;

Hweight(1,k) = interp1(fH2, H2(:,2),f(k - DCbins));

% END Hweights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% BEGIN: FluctuationStrengthBody

Hann_tmp = hanning(N, 'periodic'); 
Hann = Hann_tmp;
Hann(N/2+1:end) = 0; 

if bDebug
    % See Sontacchi1999, Abbildung 5.2 and 5.3
    figure(1)
    subplot(2,1,1); 
    plot(   1:length(Hann), Hann); 
    xlabel('Samples')
    grid on
        
    subplot(2,1,2);
    Nwin = N*2;
    hHann = freqz(Hann,1,Nwin);
    
    hBlack = freqz(Window,1,Nwin); % not used at all..
    
    f_temp = (1:Nwin)/Nwin * Fs/2;
    plot(   f_temp,20*log10(abs(hHann)) ); % figure; plot(f_temp,20*log10(abs(hBlack)));
    
    xlim([0 30]);
    xlabel('Frequency [Hz]')
    grid on
end

% Calibration between wav-level and loudness-level (assuming blackman window and FFT will follow)

Chno	=	47;
Cal	 	=	0.25;
N2		=	N/2;
q		=	1:1:N;
qb		=	N0:1:Ntop;      % N0 = 5; Ntop = 3716. N0, Ntop 1:5? % Increase N to double
freqs	=	(qb+1)*Fs/N;    % fmin = 32; ftop = 20 kHz
hBPi	=	zeros(Chno,N);
hBPrms	=	zeros(1,Chno);
mdept	=	zeros(1,Chno);
ki		=	zeros(1,Chno-2);
fi		=	zeros(1,Chno); % name changed from ri to li % Added by AO

hBPi_red = zeros(1,ceil(2*N/Fs * Fs_red)); % Added by AO
eim = zeros(Chno,2); % Added by AO

% Calculate Excitation Patterns
TempIn =  dataIn*AmpCal;
[rt,ct]=size(TempIn);
[r,c]=size(a0);

if rt~=r; 
    TempIn=TempIn';  % just transposing
    Hann = Hann'; % Added by AO 
end

TempIn	=	a0.*fft(TempIn);
Lg		=	abs(TempIn(qb));
LdB		=	To_dB(Lg);
whichL	=	find(LdB>MinExcdB);
sizL	=	length(whichL);

%% Model of Terhardt: Steepness of slopes
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

ExcAmp = zeros(sizL,47);
Slopes = zeros(sizL,47);

% For every Level above hearing threshold:
%       The masking slopes are determined along the whole frequency range (47 values)
for k=1:1:sizL
    Ltmp = LdB(whichL(k));
    Btmp = Barkno(whichL(k)+N01);

    for l = 1:1:whichZ(1,k)
        Stemp = (S1*(Btmp-(l*0.5)))+Ltmp;
        if Stemp>MinBf(l)
            Slopes(k,l)=From_dB(Stemp);
        end
    end
    for l = whichZ(2,k):1:47
        Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp;
        if Stemp>MinBf(l)
            Slopes(k,l)=From_dB(Stemp);
        end
    end
end

for k=1:1:47
    etmp = zeros(1,N);
    for l=1:1:sizL
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
    
    if bDebug
        if k == 17
            etmp_freq = etmp;
            figure(2)
            plot(freqs,20*log10(abs(etmp(qb)))), hold on
            % xlim([900 1100]), grid on
            xlabel('Frequency [Hz]')
        end
    end

    ei(k,:)	= N*real(ifft(etmp));   % This is envelope
                                    % This is envelope 
	if k == 17
        if bDebug 
            figure(3), plot(ei(k,:)), hold on
        end
    end
    
    etmp	= abs(ei(k,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is different from the Roughness calculation:

    % % DC component calculation for R:
    % h0(k)	= mean(etmp); 

    % DC component calculation for FS:
    eim(k,1) = 1/(N/4) * sum(etmp(1:N2-1));
    eim(k,2) = 1/(N/4) * sum(etmp(N2:end));

    eil = zeros(size(etmp));

    eil(1:2:end) = eim(k,2);
    eil(2:2:end) = eim(k,1);

    h0(k) = mean(abs( eil )); % Added by AO
    
    % eil = eil-h0(k); 
    eil = etmp - h0(k);
    
    % etemp - excitation patterns in frequency domain
    
    % Convolution in time before FFT
    Fei(k,:) = fft(eil ,N);
    hBPi(k,:) = 2*real(ifft(Fei(k,:))); % Added by AO % WHY 2*real?

    % % Convolution in time after synthesis, problem with subindices:
    % Fei(k,:)	= fft( eil );
    % hBPi(k,:) = conv( 2*real(ifft(Fei(k,:))),Hann); % Added by AO

    exp1 = conv(hBPi(k,:),Hann(1:N/2)); % low pass filtering, preceding downsampling
    
    % hBPi_red(k,:) = resample(exp1(1:N/2),Fs_red,Fs); % Added by AO
    tmp = resample(exp1,Fs_red,Fs); % Added by AO
    hBPi_red(k,1:length(tmp)) = tmp; % RENAME to LP
    
    % hBPi_red(k,:) = 2*real(ifft(Fei(k,:).*Hweight(1,:)));
    
    % K = N*2;
    % f0 = (1:K)/K * Fs/2;
    % figure(13);
    % plot( f0, 20*log10(abs(freqz( exp1,1,K ))) )
    % xlim([0 30])
    
    % % attempt on 10/11
    % tmp = resample( hBPi(k,:),Fs_red,Fs);
    % hBPi_red(k,1:length(tmp)) = tmp; % Added by AO
    
    opts.fs = Fs_red;
    [tmpY,tmpYdB,tmpF] = freqfft(hBPi_red(k,:)',N/2,opts);
    
    idx = 1:length(tmpF);
    tmpYf = abs(tmpY) .* Hweight(idx)'; % H - 1 x 4096 % BAND-PASS
    tmpYfdB = To_dB(tmpYf);
         
    yyyy    = 2*real(ifft(tmpYf,N));
    yy      = yyyy(1:length( hBPi_red(k,:) ));
            
    if k == 17
        if bDebug

            figure(5); 
            plot(tmpF,tmpYdB), hold on, grid on
            plot(tmpF,tmpYfdB,'r');
            plot(f_red,20*log10(abs(Hweight))+max(tmpYdB),'g');
            legend('no filter','BP filtered','BP')
            
            figure(6)
            plot(hBPi_red(k,:)'); hold on, grid on
            plot(yy,'r-')
            
            
            % figure(7); 
            % plot(freqs, 20*log10(abs(Fei(k,qb)))), hold on
            % 
            % opts.fs = Fs;
            % [ytmp tmpYdB,ftmp] = freqfft(hBPi(k,:)',N/2,opts);
            % 
            % plot(ftmp,tmpYdB,'r')
            % grid on
            % xlabel('Frequency [Hz]')
            % ylabel('Level [dB]')
            
        end
    end
    
    % hBPrms(k)       = dw_rms(hBPi_red(k,:));
    hBPrms(k)	= dw_rms(yy); % average of the band-pass filtered signals
    
    if h0(k)>0
        mdept(k) = hBPrms(k)/(h0(k));
        if mdept(k)>1
            mdept(k)=1;
        end
    else
        mdept(k)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hBPi_red(k,:) = yy'; % TEMPORAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDebug
    
    idx2plot = 14:22;
    L = size(hBPi_red,2);
    t = repmat( (1:L)/Fs_red ,length(idx2plot),1);
    BandNumber = repmat( ((idx2plot)'), 1, L);

    figure(7); 
    % mesh( t,BandNumber,hBPi(idx2plot,:) )
    mesh( t,BandNumber,abs(hBPi_red(idx2plot,:)) )
    xlabel('Time [s]')
    ylabel('Band number')
    title('Band-pass signals')
    % zlim([0 4000])
    
    idx2plot = 14:22;
    t = repmat( (1:N)/Fs ,length(idx2plot),1);
    BandNumber = repmat( ((idx2plot)'), 1, N);

    % figure(4);
    % mesh( t,BandNumber,20*log10(abs(Fei(idx2plot,:))) )
    % xlabel('Time [s]')
    % ylabel('Band number')
    % title('Band-pass signals')
    % zlim([0 4000])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find cross-correlation coefficients
for k=1:1:45
    
    cfac	=	cov(hBPi_red(k,:),hBPi_red(k+2,:));
    den	=	diag(cfac);
    den	=	sqrt(den*den');
  
    if den(2,1)>0
        ki(k)	=	cfac(2,1)/den(2,1);
    else
        ki(k)	=	0;
    end
  
end

% Calculate specific fluctuation strength fi and total fluctuation strength FS
fi(1)	=	(ki(1))^kg * (mdept(1))^p * (h0(1))^qg;
fi(2)	=	(ki(2))^kg * (mdept(2))^p * (h0(2))^qg;
for k = 3:1:45
  fi(k)	=	(ki(k-2)*ki(k))^kg * (mdept(k))^p * ( h0(k) )^qg;
end
fi(46)	=	(ki(44))^kg * (mdept(46))^p * (h0(44))^qg;
fi(47)	=	(ki(45))^kg * (mdept(47))^p * (h0(45))^qg;

FS		=	Cal*( sum(fi.^s) )^(1/s);

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
end
