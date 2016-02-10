function dataOut = FluctuationStrength_Sontacchi_vOld(insig, Fs, bDebug)
% function dataOut = FluctuationStrength_Sontacchi_vOld(insig, Fs, bDebug)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Original name: fs_offline.m
% Created on    : 26/01/2015
% Last update on: 26/01/2015 
% Last use on   : 27/01/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('Status of this code: not running')
insig = From_dB(-10)*insig;

fprintf('Average level of input signal: %.2f dB\n',rmsdb(insig)+100);

if nargin < 4
    bDebug = 0;
end

N       = 8192;
N_hop   = 4410;

Chno    = 47;

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
% dFs_red = Fs_red/N;

% Make list with Barknumber of each frequency bin
Barkno      = zeros(1,N2);
f           = N0:1:Ntop;
Barkno(f)   = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

% f_red = (1:N)/N * Fs_red;

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
ei    = zeros(Chno,N);
Fei   = zeros(Chno,N);

h0     = zeros(1,Chno);
k      = 1:1:Chno;
    
% calculate a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

method = 7; % arbitrarily, see Sontacchi1999, pp 76 of 115
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
    case 99
        kg = 0.8;   s = 1;      p = 1;      qg = 0;
end

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% H3      = interp1(fH2, H2(:,2),[1:2:60]);
% idx     = isnan(H3);
% H3(idx) = 0;
% 
% H4      = interp1(fH2, H2(:,2),[1:1:30]);
% idx     = isnan(H4);
% H4(idx) = 0;

% Hweight	= zeros(1,N);
% Hweight(1,k) = interp1(fH2, H2(:,2),f(k - DCbins));

Htmp2       = Get_Hweight_fluctuation_mod(N,Fs_red*2);
Htmp        = Get_Hweight_fluctuation_mod(N,Fs);
Hweight     = repmat( Htmp(1,:),Chno,1);
Hweight2    = repmat( Htmp2(1,:),Chno,1);

% END Hweights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Window      = blackman(N, 'periodic') .* 1.8119;
dBcorr      = 80+2.72; % originally = 80, last change to this calibration value on 26/11/2014
insig_buf   = buffer(insig, N, N-N_hop,'nodelay');
m_blocks    = size(insig_buf,2);
fi          = zeros(m_blocks,Chno);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: FluctuationStrengthBody
for idx_j = 1:m_blocks
    
    dataIn = insig_buf(:,idx_j);
    
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 

    Hann_tmp = hanning(N, 'symmetric'); 
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

        eim1(k,idx_j) = eim(k,1);
        eim2(k,idx_j) = eim(k,2);
    
    end
end

for k = 1:1:Chno
    
        sizeEim1 = size(eim1);
        eil = zeros(sizeEim1(1),sizeEim1(2)*2);
        % eil = zeros(1,sizeEim1(2)*2); % easier to debug
        eil(k,1:2:end) = eim1(k,:);
        eil(k,2:2:end) = eim2(k,:);
 
        if k == 17
            disp('');
        end
        
        h0(k) = mean(abs( eil(k,:) )); % DC value
                                        
        % h0(k) = mean(abs( eim1(k,:) )); % DC value
        %                                 % eim2 has no big influence on h0
        eil_no_DC = eil - h0(k);

        % etemp - excitation patterns in frequency domain
 
        % Convolution in time before FFT
        Fei(k,:) = fft(eil_no_DC(k,:) ,N);
        hBPi(k,:)	= 2*real( ifft(Fei(k,:).*Hweight2(k,:),N )); % 2*real(  ifft( Fei(k,:).*Hweight(k,:) ,N)  );
        h_no_filt(k,:) = 2*real( ifft(Fei(k,:)       ,N ));
        
        if k == 17
            disp('');
        end
 
        % hBPrms(k)       = rms(hBPi_red(k,:),'dim',2);
        % hBPrms(k)	= rms( eil ,'dim',2); 
        hBPrms(k)	= rms( hBPi(k,1:length(eil(k,:))) ,'dim',2); % average of the band-pass filtered signals
        if k == 17
            % rms( hBPi(k,1:length(eil(k,:))) ,'dim',2)
            % rms( h_no_filt(k,1:length(eil(k,:))) ,'dim',2)
        end
        
        if h0(k)>0
            mdept(k) = hBPrms(k)/(h0(k));
            if mdept(k)>1
                mdept(k)=1;
            end
        else
            mdept(k)=0;
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find cross-correlation coefficients
    for k=1:1:45

        cfac	=	cov(hBPi(k,1:40),hBPi(k+2,1:40));
        % cfac	=	cov(eil_no_DC(k,:),eil_no_DC(k+2,:));
        den	=	diag(cfac);
        den	=	sqrt(den*den');

        if den(2,1)>0
            ki(k)	=	abs( cfac(2,1)/den(2,1) );
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

    FS      =	Cal*( sum(fi.^s) )^(1/s);

end

% Create a cell array to return
dataOut{1} = eim1; % eim(:,1);% FS;
dataOut{2} = eim2; % eim(:,2);; %fi;
% dataOut{3} = SPL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
