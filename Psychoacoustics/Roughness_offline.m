function [R dataOut out] = Roughness_offline(insig, Fs, N, options, CParams,bDebug)
% function [R dataOut out] = Roughness_offline(insig, Fs, N, options, CParams,bDebug)
%
% 1. Description:
%       Frame-based, off-line implementation of the roughness algorithm.
% 
%       Changes by AO:
%           amp2db replaced by To_dB
%           db2amp replaced by From_dB
%           private rms renamed to dw_rms
%       
%       Outputs:
%           out - has the same format as in the PsySound toolbox
% 
% author : Matt Flax <flatmax @ http://www.flatmax.org> : Matt Flax is flatmax
% March 2006 : For the psySoundPro project
%
% revised : Farhan Rizwi, July 2007
%           Reformatted, copied and vectorised code from InitAll
%           and Hweights into the function space below.  This
%           allows us to use nested functions effeciently.
%
% contact for the original source code :
%           http://home.tm.tue.nl/dhermes/
%
% 2. Stand-alone example:
%       2.1 Unix-based example:
%           [insig fs] = Wavread('~/Documenten/MATLAB/outputs/tmp-cal/ref_rough.wav'); % Unix-based
%           [out outPsy] = Roughness_offline(insig,fs,8192);
% 
%       2.2 Multi-platform example, requires the ref_rough.wav file in the 
%           appropriate folder:
%           [insig fs] = Wavread([Get_TUe_paths('outputs') 'ref_rough.wav']); 
%           [out outPsy] = Roughness_offline(insig,fs,8192);
% 
%       2.3 Impulse response:
%           N = 8192;
%           insig  = [zeros(N/2-1,1); 1; zeros(N/2,1)];
%           fs = 44100;
%           bDebug = 1;
%           [out outPsy] = Roughness_offline(insig,fs,N,[],[],bDebug);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 10/11/2014
% Last update on: 17/01/2015 % Update this date manually
% Last use on   : 27/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    bDebug = 0;
end

if nargin < 5
    CParams = [];
    CParams.HopSize = 4096;
else
    CParams = Ensure_field(CParams,'HopSize',4096);
end

if nargin < 4
    options = [];
end

options = ef(options,'nSkipStart',0);

nSkipStart = options.nSkipStart; % to correct time series, in case an analysis frame has been excluded

N_hop   = CParams.HopSize;
insig   = insig( nSkipStart*N_hop+1:end );

if nargin < 3
    N = 8192;
end

if ~(Fs == 44100 | Fs == 40960 | Fs == 48000)
  error(['Incorrect sample rate for this roughness algorithm. Please ' ...
         're-sample original file to be Fs=44100,40960 or 48000 ' ...
         'Hz']);
end

%%%%%%%%%%%%%%%%%
% BEGIN InitAll %
%%%%%%%%%%%%%%%%%

Bark    = Get_psyparams('Bark');
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
zi      = 0.5:0.5:23.5;
zb      = sort([Bf(1,:),Cf(1,:)]);
Chno    = 47;
ei      = zeros(Chno,N);
Fei     = zeros(Chno,N);

gr = Get_psyparams('gr');

gzi    = zeros(1,Chno);
h0     = zeros(1,Chno);
k      = 1:1:Chno;
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

insig_buf   =   buffer(insig, N, N-N_hop,'nodelay');
m_blocks    =   size(insig_buf,2);
ri          =   zeros(m_blocks,Chno);

Window      = blackman(N, 'periodic') .* 1.8119;
dBcorr      = 80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_j = 1:m_blocks

    tn(idx_j) = (idx_j-1 + nSkipStart)*N_hop; % sample number to determine time
    dataIn = insig_buf(:,idx_j);
    
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); % cal to get magnitude spectrum to the power spectrum L
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calibration between wav-level and loudness-level (assuming
    % blackman window and FFT will follow)

    Cal	 	=	0.25;
    N2		=	N/2;
    q		=	1:1:N;
    qb		=	N0:1:Ntop;
    hBPi	=	zeros(Chno,N);
    hBPrms	=	zeros(1,Chno);
    mdept	=	zeros(1,Chno);
    ki		=	zeros(1,Chno-2);
    
    TempIn  = dataIn*AmpCal;
    [rt,ct] = size(TempIn);
    [r,c]   = size(a0);
    if rt~=r; TempIn=TempIn'; end

    FreqIn	= a0.*fft(TempIn);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Critical-band filterbank - Terhardt:
    [ei, etmp_fd, etmpExc_fd] = Terhardt_filterbank(FreqIn,Fs,N,qb,MinExcdB,Barkno,N01,zb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ei = zeros(  )
    for k=1:1:47
        % etmp_fd  - excitation pattern in frequency domain
        % etmp  - excitation patterns in time domain (after L242)
        % ei    - excitation patterns in time domain
        % Fei   - envelope in frequency domain, as function of fmod
        % h0    - DC component, average of half-wave rectified signal

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

    end

    if bDebug

        % params for figures 4 and 5
        idx2plot = 1:47; 

        k = idx2plot;
        
        %-- Figure 2: -----------------------------------------------------
        figure(2); 
        subplot(2,1,1)
        plot(hz2bark(freqs),20*log10(abs(etmpExc_fd(k,qb))));
        hold on, grid on
        xlabel('Critical-band rate [Bark]')
        ylabel('Excitation pattern')
        title(sprintf('Excitation pattern (Terhardt''s model), band number %.0f-%.0f',min(k),max(k)))
        ha = gca;

        subplot(2,1,2)
        plot(hz2bark(freqs),20*log10(abs(etmp_fd(k,qb)))) 
        hold on; grid on
        %xlabel('Frequency [Hz]')
        xlabel('Critical-band rate [Bark]')
        ylabel('Level [dB]')
        title(sprintf('Critical-band filterbank (Terhardt''s model), band number %.0f-%.0f',min(k),max(k)))
        %xlim([900 1100])
        ha(end+1) = gca;
        linkaxes(ha,'x')
        
        %-- Figure 4: -----------------------------------------------------
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

        %-- Figure 5: -----------------------------------------------------
        figure(5); 
        mesh( t,BandNumber,hBPi(idx2plot,:) )
        xlabel('Time [s]')
        ylabel('Band number')
        title('Band-pass signals')
        zlim([-4000 4000])

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
    ri(idx_j,1)     =	(gzi(1)*mdept(1)*ki(1))^2;
    ri(idx_j,2)     =	(gzi(2)*mdept(2)*ki(2))^2;
    for k           = 3:1:45
        ri(idx_j,k)	=	(gzi(k)*mdept(k)*ki(k-2)*ki(k))^2;
    end
    ri(idx_j,46)	=	( gzi(46)*mdept(46)*ki(44) )^2;
    ri(idx_j,47)	=	( gzi(47)*mdept(47)*ki(45) )^2;
    R(idx_j)        =	Cal*sum(ri(idx_j,:));

    SPL(idx_j) = mean(rms(dataIn));
    if SPL(idx_j) > 0
        SPL(idx_j) = To_dB(SPL(idx_j))+dBcorr+3; % -20 dBFS <--> 60 dB SPL
    else
        SPL(idx_j) = -400;
    end
    
    if bDebug
        figure; 
        plot(ri(1,:),'r'), hold on
        plot(mdept,'b-')
        plot(ki,'k')
        plot(gzi,'g--')
        grid on
        
        legend('r_i','mdepth','cc','gzi')
    end
    
end

% Create a cell array to return
dataOut{1} = R;
dataOut{2} = ri;
dataOut{3} = SPL;

out.t       = transpose(tn/Fs);
out.z       = transpose(zi);

nParam      = 1;
out.Data1   = transpose(R);
output.name{nParam} = 'Roughness';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

nParam      = 2;
out.Data2   = ri;
output.name{nParam} = 'Specific roughness';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

output.nAnalyser = 15;

out.stats.rough_tot = mean( out.Data1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% END: RoughBody

end
