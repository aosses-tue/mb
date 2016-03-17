function [R dataOut out] = Roughness_offline_debug(insig, fs, N, optsDebug)
% function [R dataOut out] = Roughness_offline_debug(insig, fs, N, optsDebug)
%
% 1. Description:
%       Off-line implementation of the roughness algorithm. It implements 
%       the algorithm using only one N-length frame.
%
% Modified by:  Alejandro Osses,
%               Matt Flax <flatmax @ http://www.flatmax.org> in March 2006 (psySoundPro project)
%
% Contact for the original source code:
%       http://home.tm.tue.nl/dhermes/
%
% 2. Stand-alone example:
%           [insig fs] = Wavread([Get_TUe_paths('outputs') delim 'ref_rough.wav']); 
%           [R outPsy out] = Roughness_offline_debug(insig(1:8192),fs,8192);
%           figure; plot(out.t,R); grid on;
%           xlabel('Time [s]'); ylabel('Roughness [asper]');
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 10/11/2014
% Last update on: 21/11/2014 
% Last use on   : 14/07/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    optsDebug = [];
end

hFig = [];

optsDebug   = ef(optsDebug,'all',0);
optsDebug   = ef(optsDebug,'ki',0);

bDebug      = optsDebug.all;
debug_ki    = optsDebug.ki;

if nargin < 3
    N = 8192;
end

if ~(fs == 44100 | fs == 40960 | fs == 48000)
  error(['Incorrect sample rate for this roughness algorithm. Please ' ...
         're-sample original file to be Fs=44100,40960 or 48000 ' ...
         'Hz']);
end

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%

Bark    = Get_psyparams('Bark');
Bark2   = [sort([Bark(:,2);Bark(:,3)]),sort([Bark(:,1);Bark(:,4)])];

N0      = round(20*N/fs)+1;
N01     = N0-1;
N50     = round(50*N/fs)-N0+1;
N2      = N/2+1;
Ntop	= round(20000*N/fs)+1;
Ntop2	= Ntop-N0+1;
dFs     = fs/N;

% Make list with Barknumber of each frequency bin
Barkno      = zeros(1,N2);
f           = N0:1:Ntop;
Barkno(f)   = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

% Make list of frequency bins closest to Cf's
Cf = ones(2,24);
for a=1:1:24
    Cf(1,a)=round(Bark((a+1),2)*N/fs)+1-N0;
    Cf(2,a)=Bark(a+1,2);
end
%Make list of frequency bins closest to Critical Band Border frequencies
Bf = ones(2,24);
Bf(1,1)=round(Bark(1,3)*N/fs);
for a=1:1:24
    Bf(1,a+1)=round(Bark((a+1),3)*N/fs)+1-N0;
    Bf(2,a)=Bf(1,a)-1;
end
Bf(2,25)=round(Bark((25),3)*N/fs)+1-N0;

%Make list of minimum excitation (Hearing Treshold)
HTres = Get_psyparams('HTres');

k = (N0:1:Ntop);
MinExcdB = interp1(HTres(:,1),HTres(:,2),Barkno(k));
  
% Initialize constants and variables
zi      = 0.5:0.5:23.5;
zb      = sort([Bf(1,:),Cf(1,:)]);
Chno    = 47;
MinBf   = MinExcdB(zb);
ei      = zeros(Chno,N);
Fei     = zeros(Chno,N);

gr = Get_psyparams('gr');

gzi    = zeros(1,Chno);
h0     = zeros(1,Chno);
k      = 1:1:Chno;
gzi(k) = sqrt(interp1(gr(1,:)',gr(2,:)',k/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates a0
a0tab =	Get_psyparams('a0tab');

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = From_dB(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

% BEGIN Hweights:

Hweight = Get_Hweight_roughness(N,fs);

% END Hweights 
%%%%%%%%%%%%%%%%

%% Stage 1, BEGIN: RoughBody
dBFS   = 100;
Window = blackman(N, 'periodic') .* 1.8119;
dBcorr = dBFS+2.72; % correction assuming that (80 dB = 0 dBFS, rms)

insig = insig .*Window;
AmpCal = From_dB(dBcorr-20)*2/(N*mean(blackman(N, 'periodic'))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration between wav-level and loudness-level (assuming
% blackman window and FFT will follow)
    
Cal	 	= 0.25;
N2		= N/2;
q		= 1:1:N;
qb		= N0:1:Ntop; % N0:1:Ntop;
hBPi	= zeros(Chno,N);
hBPrms	= zeros(1,Chno);
mdept	= zeros(1,Chno);
ki		= zeros(1,Chno-2);
ri		= zeros(1,Chno);

% Calculate Excitation Patterns
TempIn  = insig*AmpCal;
[rt,ct] = size(TempIn);
[r,c]   = size(a0);
if rt~=r; TempIn=TempIn'; end   % converts input TempIn to a column vector 1 x 8192
    
% From the time domain to the frequency domain:
FreqIn	=	a0.*fft(TempIn); 

if bDebug
    figure;
    plot((1:N)*fs/N,20*log10(abs( FreqIn )) )
    xlim([0 fs/2])
    xlabel('Frequency [Hz]')
    ylabel('Level [dB]')
    close 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Critical-band filterbank - Terhardt:
[etmp_td, etmp_fd, etmpExc_fd] = Terhardt_filterbank_debug(FreqIn,fs,N,qb,MinExcdB,Barkno,N01,zb,Chno,bDebug);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etmp_td = abs(etmp_td); % full-wave rectification

for k=1:1:Chno % each critical band number
    
    
    %% Stage 2
    
    % etmp_fd  - excitation pattern in frequency domain
    % etmp  - excitation patterns in time domain (after L242)
    % ei    - excitation patterns in time domain
    % Fei   - envelope in frequency domain
    % h0    - DC component, average of full-wave rectified signal
    
    h0(k)       = mean(etmp_td(k,:));
    Fei(k,:)	= fft( etmp_td(k,:)-h0(k) ); 
    
    hBPi(k,:)	= 2*real(  ifft( Fei(k,:).*Hweight(k,:) )  );
    hBPrms(k)	= rms(hBPi(k,:),'dim',2);
    if h0(k)>0
        mdept(k) = hBPrms(k)/h0(k);
    else
        mdept(k)=0;
    end
    
    if k==17
        disp('')
    end
    
    if bDebug
        
        figM = 3;
        figN = 2;
        plot_x = 1:size(etmp_td,2);
        plot_f = ( 1:size(etmp_fd,2) )/N*fs;
        
        if k == 15
            
            hFig(end+1) = figure(4);
            subplot(figM,figN,1)
            hp(1) = plot(plot_f,20*log10( etmp_fd(k,:) ));
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
            xlim([150/N*fs 250/N*fs]) % plotted from bin 150 to 250
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

kp(1:2)         = ki(1:2).^2;
kp(3:Chno-2)    = ki(3:Chno-2).*ki(1:Chno-4);
kp(Chno-1:Chno) = ki(Chno-3:Chno-2).^2;

if debug_ki
    
    figure;
    plot(zi,kp,'LineWidth',2); hold on
    plot(zi,mdept,'r'); grid on
    xlim([0 24]);
    
    ha = gca;
    zil = transpose( str2num( get(ha,'XTickLabel') ) );
    fil = Round(audtofreq(zil,'bark'),0);
    grid on
    set(ha,'XTickLabel',fil );
    
    xlabel('Frequency [Hz]')
    legend('kk^2','mdept')
    
end

% Calculate specific roughness ri and total roughness R
mp = min(mdept,ones(size(mdept)));

ri = (gzi.*mp.*kp).^2;
R  = Cal*sum(ri);

SPL = mean(rms(insig));
if SPL > 0
    SPL = To_dB(SPL)+dBFS; % -20 dBFS <--> 70 dB SPL
else
    SPL = -400;
end

% Create a cell array to return
dataOut{1} = R;
dataOut{2} = ri;
dataOut{3} = SPL;

out.gzi = gzi;
out.mp  = mp;
out.kp  = kp;

if bDebug
    dataOut{4} = hFig;
    disp('Total Roughness: %.2f [aspers]')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: RoughBody

end
