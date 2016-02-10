function [dataOut out] = FluctuationStrength_offline_debug(insig, Fs, N, options, bDebug)
% function [dataOut out] = FluctuationStrength_offline_debug(insig, Fs, N, options, bDebug)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
%       Comments on this implementation: cross-correlation seems not to be 
%       working properly
%
% 2. Stand-alone example:
%       r20141126_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 10/11/2014
% Last update on: 18/11/2014
% Last use on   : 29/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    bDebug = 0;
end

if nargin < 4
    options = [];
end

options = ef(options,'nSkipStart',0);

nSkipStart = options.nSkipStart; % to correct time series, in case an analysis frame has been excluded
nSkipEnd   = options.nSkipEnd;

if nargin < 3
    N = 8192*4;
end

if ~(Fs == 44100 | Fs == 40960 | Fs == 48000)
  warning(['Incorrect sample rate for this roughness algorithm. Please ' ...
           're-sample original file to be Fs=44100,40960 or 48000 Hz']);
end

N_hop   =   N/2; % 2048;
warning('temporal hop size')
insig   = insig( nSkipStart*N_hop+1:end );


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
zi      = 0.5:0.5:23.5;
zb      = sort([Bf(1,:),Cf(1,:)]);

idxtmp  = find(zb==0 | zb==1);
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

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

% BEGIN Hweights:

Htmp        = Get_Hweight_fluctuation_mod(N,Fs);
Hweight     = repmat( Htmp(1,:),Chno,1);

Hhann       = Hanning_half(8192);
% END Hweights 
%%%%%%%%%%%%%%%%

%% Stage 1, BEGIN: FluctBody

Window      = blackman(N, 'periodic') .* 1.8119;
dBcorr      = 80+2.72; % originally = 80, last change of this value on 26/11/2014

insig_buf   = buffer(insig, N, N-N_hop,'nodelay');
m_blocks    = size(insig_buf,2);
if nSkipEnd < m_blocks
    m_blocks= m_blocks - nSkipEnd;
else
    warning('nSkipEnd is larger than the amount of analysis frames...')
end
fi          =   zeros(m_blocks,Chno);

if bDebug 
    if m_blocks > 1
        warning('Only first frame is going to be plotted')
        pause(2)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_j = 1:m_blocks

    tn(idx_j) = (idx_j-1 + nSkipStart)*N_hop; % sample number to determine time
    dataIn = insig_buf(:,idx_j);
    
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 
    
    % Calibration between wav-level and loudness-level (assuming
    % blackman window and FFT will follow)

    Cal	 	=	0.25; % 0.138;
    N2		=	N/2;
    q		=	1:1:N;
    qb		=	N0:1:Ntop;
    freqs	=	(qb+1)*Fs/N;
    hBPi	=	zeros(Chno,N);
    hCrossi	=	zeros(Chno,N);
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

        if k == 17
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
        ei(k,:)     = N*real(ifft(etmp));
        etmp_td(k,:)= abs(ei(k,:));
        
        h0(k)       = mean(etmp_td(k,:));
        Fei(k,:)	= fft( etmp_td(k,:)-h0(k) ); % changes the phase but not the amplitude
        
        hBPi(k,:)	= 2*real(  ifft( Fei(k,:).*Hweight(k,:) ,N)  );
        hCrossi(k,:)= 2*real(  ifft( Fei(k,:)               ,N)  );
        
        hBPrms(k)	= rms(hBPi(k,:),'dim',2);
        
       
        if k == 17

            if bDebug
                optstm.fs = 44100;
                [tmy tmydB1 ff1] = freqfft(ei(k,:)',N/2,optstm);

                optstm.fs = 44100;
                [tmy3 tmydB3 ff3] = freqfft(Hhann,N/2,optstm);

                figure;
                plot(ff1,tmydB1,ff3,tmydB3)
                legend('ei','win')
                xlim([0 150])
            end

        end
        
        if h0(k)>0
            
            mdept(k) = hBPrms(k)/h0(k);
            % mdeptr(k) = hBPrmsc(k)/h0c(k);
            
            if mdept(k)>1
                % mdept(k)=1;
            end
            % if mdeptr(k)>1
            %     mdeptr(k)=1;
            % end

        else
            mdept(k)=0;
            % mdeptr(k) = 0;
        end

        if bDebug

            figM = 3;
            figN = 2;
            plot_x = 1:size(etmp_td,2);
            plot_f = ( 1:size(etmp_fd,2) )/N*Fs;

            if k == 15

                hFig(end+1) = figure(4);
                subplot(figM,figN,1)
                hp(1) = plot(plot_f,abs(etmp_fd(k,:)));
                ha = gca;
                grid on, hold on
                title(sprintf('Excitation patterns in freq. domain\nNum band = %.0f',k))
                xlabel('Frequency [Hz]')

                subplot(figM,figN,2)
                hp(1) = plot(plot_x,abs(etmp_td(k,:))); grid on, hold on
                hp(2) = plot(minmax(plot_x),[h0(k) h0(k)],'r');
                hp(3) = plot(plot_x,hBPi(k,:),'k','LineWidth',2);

                legend(hp(2:3),   {sprintf( 'h0=%.1f',h0(k)),...
                                   sprintf( 'hBPrms = %.1f',hBPrms(k))});
                ha1 = gca;

                title( sprintf('Excitation patterns in time domain (half rectified)\nNum band = %.0f, mdepth = %.2f',k,mdept(k)) )
                xlabel('Samples')

            end

            if k == 17
                figure(4)
                subplot(figM,figN,3)
                plot(plot_f,abs(etmp_fd(k,:)))
                ha(end+1) = gca;
                grid on 
                title(sprintf('Num band = %.0f',k))
                xlabel('Frequency [Hz]')

                subplot(figM,figN,4)
                hp(1) = plot(plot_x,abs(etmp_td(k,:))); grid on, hold on
                hp(2) = plot(minmax(plot_x),[h0(k) h0(k)],'r');
                hp(3) = plot(plot_x,hBPi(k,:),'k','LineWidth',2);

                legend(hp(2:3),   {sprintf( 'h0=%.1f',h0(k)),...
                                   sprintf( 'hBPrms = %.1f',hBPrms(k))});
                ha1(end+1) = gca;

                title(sprintf('Num band = %.0f, mdepth = %.2f',k,mdept(k)))
                xlabel('Samples')
            end

            if k == 18
                figure(4)
                subplot(figM,figN,5)
                plot(plot_f,abs(etmp_fd(k,:)))
                ha(end+1) = gca;
                linkaxes(ha,'xy');
                ylim([0 8000])
                xlim([850 1300]) 
                grid on
                title(sprintf('Num band = %.0f',k))
                xlabel('Frequency [Hz]')

                subplot(figM,figN,6)
                hp(1) = plot(plot_x,abs(etmp_td(k,:))); grid on, hold on
                hp(2) = plot(minmax(plot_x),[h0(k) h0(k)],'r');
                hp(3) = plot(plot_x,hBPi(k,:),'k','LineWidth',2);

                legend(hp(2:3),   {sprintf( 'h0=%.1f',h0(k)),...
                                   sprintf( 'hBPrms = %.1f',hBPrms(k))});
                ha1(end+1) = gca;

                linkaxes(ha1,'xy')
                ylim([-2000 40000])
                xlim(minmax(plot_x))
                grid on, hold on
                plot(minmax(plot_x),[h0(k) h0(k)],'r')
                plot(plot_x,hBPi(k,:),'k','LineWidth',2)
                title(sprintf('Num band = %.0f, mdepth = %.2f',k,mdept(k)))
                xlabel('Samples')

                disp('')
            end

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

    for k=1:1:Chno-2
        ccfac	=	cov(hCrossi(k,:),hCrossi(k+2,:));
        denc	=	diag(ccfac);
        denc	=	sqrt(denc*denc');
        if denc(2,1)>0
            kki(k)	=	ccfac(2,1)/denc(2,1);
        else
            kki(k)	=	0;
        end
    end
    
    % Calculate specific fluctuation strength fi and total FS
    fi(idx_j,1)         =	(h0(1)^qg)*(mdept(1)^p)*(ki(1)^kg);
    fi(idx_j,2)         =	(h0(2)^qg)*(mdept(2)^p)*(ki(2)^kg);

    for k = 3:1:Chno-2
        fi(idx_j,k)     =	(h0(k)^qg)*(mdept(k)^p)*( ki(k-2)*ki(k) )^kg;
    end

    fi(idx_j,Chno-1)	=	(h0(Chno-1)^qg)*(mdept(Chno-1)^p)*(ki(Chno-3)^kg); % ki
    fi(idx_j,Chno  )	=	(h0(Chno  )^qg)*(mdept(Chno  )^p)*(ki(Chno-2)^kg);
    FS(idx_j)           =	Cal*( sum(fi(idx_j,:)) )^(1/s);

    SPL(idx_j) = mean(rms(dataIn));
    if SPL(idx_j) > 0
        SPL(idx_j) = To_dB(SPL(idx_j))+90; 
    else
        SPL(idx_j) = -400;
    end
    
    if bDebug
        figure; 
        plot(fi(1,:),'r'), hold on
        plot(mdept,'b-')
        plot(ki,'k')
        grid on
        title(sprintf('FS = %.2f [vacil]',FS(1)))
        legend('f_i','mdepth','cc')
    end
    
    bDebug = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: FluctBody

% Create a cell array to return
dataOut{1} = FS;
dataOut{2} = fi;
dataOut{3} = SPL;

out.t       = transpose(tn/Fs);
out.z       = transpose(zi);

nParam      = 1;
out.Data1   = transpose(FS);
output.name{nParam} = 'Fluctuation strength';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

nParam      = 2;
out.Data2   = fi;
output.name{nParam} = 'Specific fluctuation strength';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

output.nAnalyser = 20;

end
