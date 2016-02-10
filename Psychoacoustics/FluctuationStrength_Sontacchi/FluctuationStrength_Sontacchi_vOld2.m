function dataOut = FluctuationStrength_Sontacchi_vOld2(insig, Fs, bDebug)
% function dataOut = FluctuationStrength_Sontacchi_vOld2(insig, Fs, bDebug)
%
% 1. Description:
%       Off-line implementation of the Fluctuation Strength algorithm based 
%       on Roughness model.
% 
%       Changes introduced by AO:
%           - amp2db replaced by To_dB
%           - db2amp replaced by From_dB
%           - private rms renamed to dw_rms
%           - Compatible with any Fs. If you use Fs smaller than 44.100 Hz, 
%             then the calculation will assume no contribution for all the µ
%             bands above Fs/2. 
%
% author : Matt Flax <flatmax @ http://www.flatmax.org> : Matt Flax is flatmax
% March 2006 : For the psySoundPro project
%
% revised : Farhan Rizwi
%           July '07
%           Reformatted, copied and vectorised code from InitAll and Hweights 
%           into the function space below. This allows us to use nested 
%           functions effeciently.
% contact for the original source code :
% http://home.tm.tue.nl/dhermes/
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Original file name: FluctuationStrength_debug.m
% Created on    : 10/11/2014
% Last update on: 18/11/2014 
% Last use on   : 18/11/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    bDebug = 0;
end

N_insig = length(insig);
N = 8192;

if ~(Fs == 44100 | Fs == 40960 | Fs == 48000)
  warning(['Incorrect sample rate for this roughness algorithm. Please ' ...
           're-sample original file to be Fs=44100,40960 or 48000 Hz']);
end

%% BEGIN InitAll %

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

h0     = zeros(1,Chno);
k      = 1:1:Chno;
    
% calculate a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 9; % similar setting than Roughness, see Sontacchi1998, pp 76 of 115
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

hFig = [];
kref = 25;

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

% BEGIN Hweights:

N_hop       =   2048*2;
Htmp        = Get_Hweight_fluctuation_mod(N,Fs/N_hop); % correct other half
Hweight     = repmat( Htmp(1,:),Chno,1);

Hhann       = Hanning_half(8192);
% END Hweights 
%%%%%%%%%%%%%%%%

%% Stage 1, BEGIN: FluctBody

Window      =   blackman(N, 'periodic') .* 1.8119;
dBcorr      =   80+10.72; % originally = 80

insig_buf   =   buffer(insig, N, N-N_hop,'nodelay');
m_blocks    =   size(insig_buf,2);
fi          =   zeros(1,Chno);
eim         =   zeros(Chno,m_blocks*2);
Fei         =   zeros(Chno,N);

if bDebug 
    if m_blocks > 1
        warning('Only first frame is going to be plotted')
        pause(2)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_j = 1:m_blocks

    fprintf('Processing block %.0f out of %.0f\n',idx_j,m_blocks);
    
    dataIn = insig_buf(:,idx_j);
    
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 
    
    % Calibration between wav-level and loudness-level (assuming
    % blackman window and FFT will follow)

    Cal	 	=	0.138; % 0.25;
    N2		=	N/2;
    q		=	1:1:N;
    qb		=	N0:1:Ntop;
    freqs	=	(qb+1)*Fs/N;
    hBPi	=	zeros(Chno,m_blocks*2);
    hBPrms	=	zeros(1,Chno);
    mdept	=	zeros(1,Chno);
    ki		=	zeros(1,Chno-2);

    % Calculate Excitation Patterns
    TempIn  = dataIn*AmpCal;    % Just windowed temporal signal
    [rt,ct] = size(TempIn);
    [r,c]   = size(a0);
    if rt  ~= r; TempIn=TempIn'; end

    TempIn	= a0.*fft(TempIn);
    Lg		= abs(TempIn(qb));
    LdB		= To_dB(Lg);
    whichL	= find(LdB>MinExcdB); % Frequency components above threshold
    sizL	= length(whichL);

    if bDebug
    
        hFig = figure(2);

        subplot(2,1,1)
        semilogx(freqs,MinExcdB), hold on
        plot(freqs(whichL),LdB(whichL),'r')

        xlabel('Frequency [Hz]')
        ylabel('Magnitude')

        legend('Hearing threshold','components above threshold')

    end
    
    % Assessment of slopes (Terhardt)
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
    % End: assessment of slopes (Terhardt)
    
    if bDebug
    
        figure(2);
        subplot(2,1,2)

        idxq = whichL;
        L = length(idxq);

        idx2plot = 14:22;
        BandNumber = repmat( ((idx2plot)'),1,L);
        freqsm = repmat(freqs(whichL),length(idx2plot),1);
        mesh(freqsm',BandNumber',Slopes(:,idx2plot)), hold on

        xlabel('Frequency [Hz]')
        ylabel('CB number')
        title('Determined masking patterns')

    end

    for k=1:1:Chno

        etmp = zeros(1,N);

        for l=1:1:sizL % each level above threshold

            N1tmp = whichL(l);
            N2tmp = N1tmp + N01;
            if (whichZ(1,l) == k)
                ExcAmp(l, k) = 1;
            elseif (whichZ(2,l) == k)
                ExcAmp(l, k) = 1;
            elseif (whichZ(2,l) > k)
                ExcAmp(l,k) = Slopes(l,k+1)/Lg(N1tmp);
            else
                ExcAmp(l,k) = Slopes(l,k-1)/Lg(N1tmp);
            end
            etmp(N2tmp) = ExcAmp(l,k)*TempIn(N2tmp); 
        end

        if k == kref
            if bDebug
                hFig(end+1) = figure(3); 
                plot(freqs,20*log10(abs(etmp(qb))),'LineWidth',2), hold on
                plot(freqs,20*log10(abs(TempIn(qb))),'r--')
                hold on; grid on
                xlabel('Frequency [Hz]')
                ylabel('Level [dB]')
                title('Spectral components: Terhardt''s model')
                legend('only excitation','M + all freq components')
                xlim([900-200 1100+200])
            end
        end
        
        %% Stage 2
        
        % etmp_fd  - excitation pattern in frequency domain
        % etmp  - excitation patterns in time domain (after L242)
        % ei    - excitation patterns in time domain
        % Fei   - envelope in frequency domain, as function of fmod
        % h0    - DC component, average of half-wave rectified signal

        etmp_fd(k,:)= etmp;
        % ei(k,:)     = N*real(ifft(etmp .* fft(Hhann)'.*Hweight(k,:) ));
        ei(k,:)     = N*real(ifft(etmp .*fft(Hhann)'));
        % eiconv(k,:) = conv(ei(k,:),Hhann);
        
        % etmp_td(k,:)= abs(ei(k,:));
        
        eim_idx1 = 2*idx_j - 1;
        eim_idx2 = 2*idx_j;
        
        % eim(k,eim_idx1)     = 1/(N/4) * sum(abs(etmp(1:N2-1)));
        % eim(k,eim_idx2)     = 1/(N/4) * sum(abs(etmp(N2:end)));
        eim(k,eim_idx1)     = 1/(N/4) * sum(abs(ei(k,1:N2-1)));
        eim(k,eim_idx2)     = 1/(N/4) * sum(abs(ei(k,N2:end)));
        
        if k == kref
            disp('')
        end
        
    end
end
    
for k=1:1:Chno
        
    % DC component calculation for FS:
    h0(k)       = mean(eim(k,:));
    Fei(k,:)	= fft( eim(k,:)-h0(k) ,N); % changes the phase but not the amplitude
        
    hBPitmp     = 2*real(  ifft( Fei(k,:).*Hweight(k,:) ,N)  );
    hBPi(k,:)   = hBPitmp(1:size(eim,2));
    %hBPi(k,:)   = eim(k,:)-h0(k);
    hBPrms(k)	= dw_rms(hBPi(k,:));
         
    if k == kref
        % ktmp = 25; figure; plot(hBPi(ktmp,:)), hold on, plot(hBPi(ktmp+10,:),'r');
        disp('')
    end

    if h0(k)>0
        
        mdept(k) = hBPrms(k)/h0(k);
        if mdept(k)>1
            mdept(k)=1;
        end
        
    else
        
        mdept(k)=0;
        
    end
        
end
    
    %% Stage 3
    
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
    fi(1)         =	(h0(1)^qg)*(mdept(1)^p)*(ki(1)^kg);
    fi(2)         =	(h0(2)^qg)*(mdept(2)^p)*(ki(2)^kg);
     
    for k = 3:1:Chno-2
        fi(k)     =	(h0(k)^qg)*(mdept(k)^p)*( ki(k-2)*ki(k) )^kg;
    end
     
    fi(Chno-1)	=	(h0(Chno-1)^qg)*(mdept(Chno-1)^p)*(ki(Chno-3)^kg);
    fi(Chno  )	=	(h0(Chno  )^qg)*(mdept(Chno  )^p)*(ki(Chno-2)^kg);
    FS           =	Cal*( sum(fi) )^(1/s);
     
    SPL = mean(rms(insig));
    if SPL > 0
        SPL = To_dB(SPL)+dBcorr+3; % -20 dBFS <--> 60 dB SPL
    else
        SPL = -400;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: FluctBody

% Create a cell array to return
dataOut{1} = FS;
dataOut{2} = fi;
dataOut{3} = SPL;

end
