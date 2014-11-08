function out = fluctuation(Fs, flag)
% function out = fluctuation(Fs, flag)
%
% 1. Description:
%       Algorithm to determine Fluctuation Strength.
%       If nargin == 2, the window length is returned.
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
%
% contact for the original source code :
% http://home.tm.tue.nl/dhermes/
%
% the following files are support files for this code :
%                InitFs44100N8192.mat
%                InitFs40960N8192.mat
%                InitFs48000N8192.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 8192;
if ~(Fs == 44100 | Fs == 40960 | Fs == 48000)
  error(['Incorrect sample rate for this roughness algorithm. Please ' ...
         're-sample original file to be Fs=44100,40960 or 48000 ' ...
         'Hz']);
end
if nargin == 2
  out = N;
  return
end

% Create and return a function handle to the roughness routine
out = @FluctBody;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute the fluctuation strength algorithm
function dataOut = FluctBody(dataIn)
    
    % Added by AO:
    Hann = hanning(N, 'periodic'); 
    Hann(N/2+1:end) = 0; % subplot(2,1,1); plot(Hann); figure; freqz(Hann,1,4096);
    % end Added by AO
    
    AmpCal = db2amp(80)*2/(N*mean(blackman(N, 'periodic')));
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
    Fs_red  = 20; % in Hz, Added by AO
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
    LdB		=	amp2db(Lg);
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
                Slopes(k,l)=db2amp(Stemp);
            end
        end
        for l = whichZ(2,k):1:47
            Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp;
            if Stemp>MinBf(l)
                Slopes(k,l)=db2amp(Stemp);
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
            etmp(N2tmp) = ExcAmp(N1tmp,k)*TempIn(N2tmp); % plot(freqs,20*log10(abs(etmp(qb)))), hold on
        end
        
        if k == 17
            disp('')
        end
        ei(k,:)	= N*real(ifft(etmp));   % This is envelope
                                        % This is envelope % figure(2), plot(ei(k,:)), hold on
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
        
        eil = eil-h0(k);
        % Fei(k,:)	= fft( eil );
        % hBPi(k,:)	= 2*real(ifft(Fei(k,:).*Hweight(k,:))); % commented by AO
        
        % Convolution in time before FFT
        % Fei(k,:) = fft( conv(eil,Hann) ,N); % Be careful, already convolved!!
        Fei(k,:) = fft(eil ,N);
        hBPi(k,:) = 2*real(ifft(Fei(k,:))); % Added by AO
        
        % % Convolution in time after synthesis, problem with subindices:
        % Fei(k,:)	= fft( eil );
        % hBPi(k,:) = conv( 2*real(ifft(Fei(k,:))),Hann); % Added by AO
        
        hBPi_red(k,:) = resample(conv(hBPi(k,:),Hann),Fs_red,Fs); % Added by AO
        % L = lenght(hBPi_red); % Added by AO
        
                                    % figure(2), plot(2*real(ifft(Fei(k,:))),'r')
                                    % figure(2), plot(2*real(ifft(Fei(k,:)).*Hann),'g')
                                        
                                    % figure(2), plot(  conv( 2*real(ifft(Fei(k,:))),Hann),'g')
        
        % figure(10); freqz(hBPi(k,:),1,4096)
        % figure(11); freqz(hBPr,1,4096)
                
        %hBPrms(k)	= rms(hBPi(k,:));
        hBPrms(k)	= rms(hBPi_red(k,:));
        if h0(k)>0
            mdept(k) = hBPrms(k)/h0(k);
            if mdept(k)>1
                mdept(k)=1;
            end
        else
            mdept(k)=0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % find cross-correlation coefficients
    for k=1:1:45
      % cfac	=	cov(hBPi(k,:),hBPi(k+2,:));
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
      SPL = amp2db(SPL)+83; % -20 dBFS <--> 60 dB SPL
    else
      SPL = -400;
    end
    
    % Create a cell array to return
    dataOut{1} = FS;
    dataOut{2} = fi;
    dataOut{3} = SPL;
  end % RoughBody
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%
Bark = [0     0	   50	 0.5
        1   100	  150	 1.5
        2   200	  250	 2.5
        3   300	  350	 3.5
        4   400	  450	 4.5
        5   510	  570	 5.5
        6   630	  700	 6.5
        7   770	  840	 7.5
        8   920	 1000	 8.5
        9  1080	 1170	 9.5
        10  1270 1370	10.5
        11  1480 1600	11.5
        12  1720 1850	12.5
        13  2000 2150	13.5
        14  2320 2500	14.5
        15  2700 2900	15.5
        16  3150 3400	16.5
        17  3700 4000	17.5
        18  4400 4800	18.5
        19  5300 5800	19.5
        20  6400 7000	20.5
        21  7700 8500	21.5
        22  9500 10500	22.5
        23 12000 13500	23.5
        24 15500 20000	24.5];

Bark2	= [sort([Bark(:,2);Bark(:,3)]),sort([Bark(:,1);Bark(:,4)])];
fi = 1; % fi = 20 for roughness
N0	= round(fi*N/Fs)+1;
N01	= N0-1;
N50     = round(50*N/Fs)-N0+1;
N2	= N/2+1;
Ntop	= round(20000*N/Fs)+1;
Ntop2	= Ntop-N0+1;
dFs	= Fs/N;

% Make list with Barknumber of each frequency bin
Barkno	  = zeros(1,N2);
f	  = N0:1:Ntop;
Barkno(f) = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

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
HTres= [	0		130
            0.01   70
            0.17	 60
            0.8	 30
            1		 25
            1.5	 20
            2		 15
            3.3	 10
            4		  8.1
            5		  6.3
            6		  5
            8		  3.5
            10		  2.5
            12		  1.7
            13.3	  0
            15		 -2.5
            16		 -4
            17		 -3.7
            18		 -1.5
            19		  1.4
            20		  3.8
            21		  5
            22		  7.5
            23 	 15
            24 	 48
            24.5 	 60
            25		130];

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
a0tab =	[ 0	 0
          10	 0
          12	 1.15
          13	 2.31
          14	 3.85
          15	 5.62
          16	 6.92
          16.5	 7.38
          17	 6.92
          18	 4.23
          18.5	 2.31
          19	 0
          20	-1.43
          21	-2.59
          21.5	-3.57
          22	-5.19
        22.5	-7.41
          23	-11.3
          23.5	-20
          24	-40
          25	-130
          26	-999];

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = db2amp(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

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

end % end fluctuation strength