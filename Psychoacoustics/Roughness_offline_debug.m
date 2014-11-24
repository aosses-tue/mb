function dataOut = Roughness_offline_debug(dataIn, Fs, N, bDebug)
% function dataOut = Roughness_offline_debug(dataIn, Fs, N, bDebug)
%
% 1. Description:
%       Off-line implementation of the roughness algorithm.
%       Corrections by AO:
% 
%           ExcAmp(N1tmp, k) changed by ExcAmp(l, k) - this optimises the 
%           memory allocation reserved for ExpAmp
% 
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
% Last update on: 21/11/2014 % Update this date manually
% Last use on   : 21/11/2014 % Update this date manually
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

%% BEGIN InitAll %

Bark = Get_psyparams('Bark');
Bark2   = [sort([Bark(:,2);Bark(:,3)]),sort([Bark(:,1);Bark(:,4)])];

N0      = round(20*N/Fs)+1;
N01     = N0-1;
N50     = round(50*N/Fs)-N0+1;
N2      = N/2+1;
Ntop	= round(20000*N/Fs)+1;
Ntop2	= Ntop-N0+1;
dFs     = Fs/N;

% Make list with Barknumber of each frequency bin
Barkno      = zeros(1,N2);
f           = N0:1:Ntop;
Barkno(f)   = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

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
HTres = Get_psyparams('HTres');

k = (N0:1:Ntop);
MinExcdB = interp1(HTres(:,1),HTres(:,2),Barkno(k));
  
% Initialize constants and variables
zi    = 0.5:0.5:23.5;
zb    = sort([Bf(1,:),Cf(1,:)]);
MinBf = MinExcdB(zb);
ei    = zeros(47,N);
Fei   = zeros(47,N);

gr = Get_psyparams('gr');

gzi    = zeros(1,47);
h0     = zeros(1,47);
k      = 1:1:47;
gzi(k) = sqrt(interp1(gr(1,:)',gr(2,:)',k/2));
    
% calculate a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

% BEGIN Hweights:

Hweight = Get_Hweight_roughness(N,Fs);

% END Hweights 
%%%%%%%%%%%%%%%%

%% Stage 1, BEGIN: RoughBody

Window = blackman(N, 'periodic') .* 1.8119;
dBcorr = 80+10.72; % originally = 80

dataIn = dataIn .*Window;
AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration between wav-level and loudness-level (assuming
% blackman window and FFT will follow)
    
Chno	=	47;
Cal	 	=	0.25;
N2		=	N/2;
q		=	1:1:N;
qb		=	N0:1:Ntop; % N0:1:Ntop;
freqs	=	qb*Fs/N; % freqs	=	(qb+1)*Fs/N;
hBPi	=	zeros(Chno,N);
hBPrms	=	zeros(1,Chno);
mdept	=	zeros(1,Chno);
ki		=	zeros(1,Chno-2);
ri		=	zeros(1,Chno);

% Calculate Excitation Patterns
TempIn  =   dataIn*AmpCal;
[rt,ct] =   size(TempIn);
[r,c]   =   size(a0);
if rt~=r; TempIn=TempIn'; end   % converts input TempIn to a column vector 1 x 8192
    
% From the time domain to the frequency domain:
TempIn	=	a0.*fft(TempIn);    

if bDebug
    figure;
    plot((1:N)*Fs/N,20*log10(abs( TempIn )) )
    xlim([0 Fs/2])
    xlabel('Frequency [Hz]')
    ylabel('Level [dB]')
    close 
end

Lg		=	abs(TempIn(qb));    % It takes only the frequencies of interest
                                % semilogx(freqs,Lg), xlabel('Frequency [Hz]'), ylabel('Magnitude')
LdB		=	To_dB(Lg);
whichL	=	find(LdB>MinExcdB); % Frequency components above hearing thresholds
sizL	=	length(whichL);

if bDebug
    
    hFig = figure(2);
    
    subplot(2,1,1)
    semilogx(freqs,MinExcdB), hold on, grid on
    plot(freqs(whichL),LdB(whichL),'r')
    
    xlim([minmax(freqs(whichL))].*[0.5 2]) % One octave below and above the frequencies above Thres.
    ylim([min(MinExcdB) max(max(LdB),max(MinExcdB))*2])
    xlabel('Frequency [Hz]')
    ylabel('Magnitude')
    
    legend('Hearing threshold','components above threshold')
    
end

% Assessment of slopes (Terhardt) for freq components above thres.
S1 = -27;
S2 = zeros(1,sizL);

for w = 1:1:sizL;
    % Steepness of upper slope [dB/Bark] in accordance with Terhardt
    steep = -24-(230/freqs(w))+(0.2*LdB(whichL(w)));
    if steep < 0
        S2(w) = steep;
    end
end
% END: Assessment of slopes

whichZ      = zeros(2,sizL);
qd          = 1:1:sizL;

% In which critical band number are the levels above threshold located:
whichZ(1,:)	= floor(2*Barkno(whichL(qd)+N01));  % idxs up to which S1 is going to be used
whichZ(2,:)	= ceil(2*Barkno(whichL(qd)+N01));   % idxs from which S2 is going to be used

ExcAmp = zeros(sizL,47);
Slopes = zeros(sizL,47);

% One masking pattern (as a function of freq) per each level LdB above 
% threshold. Therefore sizL patterns will be obtained.
% The masking patterns are going to be in the critical bands between
% whichZ(1,:) and whichZ(2,:)
for k=1:1:sizL
    Ltmp = LdB(whichL(k));          % Level [dB]
    Btmp = Barkno(whichL(k)+N01);   % Frequency [Bark]

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
% End: assessment of slopes (Terhardt)

if bDebug
    
    figure(2);
    subplot(2,1,2)
    
    idxq = whichL;
    L = length(idxq);
    
    idx2plot = 14:22;
    BandNumber = repmat( ((idx2plot)'),1,L);
    
    lvl_idx = 1:length(freqs(whichL));
    
    level_sm = repmat(lvl_idx,length(idx2plot),1);
    mesh(level_sm',BandNumber',To_dB( Slopes(:,idx2plot) )), hold on
    
    xlim(minmax(lvl_idx))
    
    xlabel(sprintf('Level_{idx}\nnumber of levels above threshold = %.0f',max(lvl_idx)))
    ylabel('CB number')
    title('Determined masking patterns')
    
end

if bDebug
    etmp_array   = zeros(47,N);
end

for k=1:1:47 % each critical band number
    
    etmp    = zeros(1,N);
    
    for l=1:1:sizL % each level above threshold
        
        N1tmp = whichL(l);
        N2tmp = N1tmp + N01; % N1tmp 'corrected' to match Lg idxs with TempIn idxs
        
        if (whichZ(1,l) == k)  % If lower limit is equal to k (Second enters here)
            ExcAmp(l, k) = 1;
            if bDebug; idxExcAmp(l,k)=1; end
            
        elseif (whichZ(2,l) == k) % If upper limit is equal to k (Third enters here)
            ExcAmp(l, k) = 1;
            if bDebug; idxExcAmp(l,k)=2; end
            
        elseif (whichZ(2,l) > k) % First enters here
            ExcAmp(l,k) = Slopes(l,k+1)/Lg(N1tmp);  % Increasing slopes, increasing Lg as k increases:
                                                    % So: ExpAmp increases with values from 0 to 1
            if bDebug; idxExcAmp(l,k)=3; end
            
        else % Fourth enters here
            ExcAmp(l,k) = Slopes(l,k-1)/Lg(N1tmp);  % Decreasing slopes, decreasing Lg as k decreases
                                                    % So: ExpAmp decreases with values from 1 to 0
            if bDebug; idxExcAmp(l,k)=4; end
            
        end
        
        if bDebug; disp(sprintf('CB number k = %.0f, ExcAmp=%.2f, entered ''if'' %.0f',k,ExcAmp(l, k),idxExcAmp(l,k))); end
        
        % figure; ktmp = 13; plot(Slopes(:,ktmp)./transpose(Lg(whichL)));
        % figure; plot(freqs, abs(etmp(qb)) ), hold on; plot(freqs(whichL),ExcAmp(:,ktmp)/max(ExcAmp(:,ktmp))*max(abs(etmp)),'r'), xlim(minmax(freqs(whichL)).*[0.8 1.2]), pause(2), close
        
    end
    etmp( whichL+N01 ) = ExcAmp(:,k)'.*TempIn( whichL+N01 );
    if bDebug
        etmp_array(k, whichL+N01 ) = ExcAmp(:,k)'.*TempIn( whichL+N01 ); 
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
            legend('Excitation, etmp (to be used for calculation)','M + all freq components')
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
    
    hBPi(k,:)	= 2*real(  ifft( Fei(k,:).*Hweight(k,:) )  );
    hBPrms(k)	= dw_rms(hBPi(k,:));
    
    if h0(k)>0
        mdept(k) = hBPrms(k)/h0(k);
        if mdept(k)>1
            mdept(k)=1;
        end
    else
        mdept(k)=0;
    end
    
    if bDebug
        
        figM = 3;
        figN = 2;
        plot_x = 1:size(etmp_td,2);
        plot_f = ( 1:size(etmp_fd,2) )/N*Fs;
        
        if k == 15
            
            hFig(end+1) = figure(4);
            subplot(figM,figN,1)
            hp(1) = plot(plot_f,20*log10(abs(etmp_fd(k,:))));
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
            plot(plot_f,20*log10(abs(etmp_fd(k,:))));
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
            plot(plot_f,20*log10(abs(etmp_fd(k,:))));
            ha(end+1) = gca;
            linkaxes(ha,'xy');
            ylim([0 2500])
            xlim([150/N*Fs 250/N*Fs]) % plotted from bin 150 to 250
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
            ylim([-2000 12000])
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
for k = 3:1:45
  ri(k)	=	(gzi(k)*mdept(k)*ki(k-2)*ki(k))^2;
end
ri(46)	=	(gzi(46)*mdept(46)*ki(44))^2;
ri(47)	=	(gzi(47)*mdept(47)*ki(45))^2;
R		=	Cal*sum(ri);

SPL = mean(rms(dataIn));
if SPL > 0
    SPL = To_dB(SPL)+dBcorr+3; % -20 dBFS <--> 60 dB SPL
else
    SPL = -400;
end

% Create a cell array to return
dataOut{1} = R;
dataOut{2} = ri;
dataOut{3} = SPL;

if bDebug
    dataOut{4} = hFig;
    disp('Total Roughness: %.2f [aspers]')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: RoughBody

end
