function [dataOut out] = Roughness_offline(insig, Fs, N, options, CParams,bDebug)
% function [dataOut out] = Roughness_offline(insig, Fs, N, options, CParams,bDebug)
%
% 1. Description:
%       Frame-based, off-line implementation of the roughness algorithm.
%       Changes:
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
%           [insig fs] = Wavread('~/Documenten/MATLAB/outputs/tmp-cal/ref_rough.wav');
%           [insig fs] = Wavread([Get_TUe_paths('outputs') 'ref_rough.wav']);
%           [out outPsy] = Roughness_offline(insig,fs,8192,0);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 10/11/2014
% Last update on: 17/01/2015 % Update this date manually
% Last use on   : 22/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    bDebug = 0;
end

if nargin < 5
    CParams = [];
    CParams.HopSize = 4096;
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
zi      = 0.5:0.5:23.5;
zb      = sort([Bf(1,:),Cf(1,:)]);
MinBf   = MinExcdB(zb);
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
dBcorr      = 81; % originally = 80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_j = 1:m_blocks

    tn(idx_j) = (idx_j-1 + nSkipStart)*N_hop; % sample number to determine time
    dataIn = insig_buf(:,idx_j);
    
    dataIn = dataIn .*Window;
    AmpCal = From_dB(dBcorr)*2/(N*mean(blackman(N, 'periodic'))); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calibration between wav-level and loudness-level (assuming
    % blackman window and FFT will follow)

    Cal	 	=	0.25;
    N2		=	N/2;
    q		=	1:1:N;
    qb		=	N0:1:Ntop;
    freqs	=	(qb+1)*Fs/N;
    hBPi	=	zeros(Chno,N);
    hBPrms	=	zeros(1,Chno);
    mdept	=	zeros(1,Chno);
    ki		=	zeros(1,Chno-2);
    
    % Calculate Excitation Patterns
    TempIn  = dataIn*AmpCal;
    [rt,ct] = size(TempIn);
    [r,c]   = size(a0);
    if rt~=r; TempIn=TempIn'; end

    TempIn	=	a0.*fft(TempIn);
    Lg		=	abs(TempIn(qb));
    LdB		=	To_dB(Lg);
    idx_L	=	find(LdB>MinExcdB); % index of levels above threshold
    sizL	=	length(idx_L);

    % steepness of slopes (Terhardt)
    S1 = -27;
    S2 = zeros(1,sizL);

    tmpF = [];
    tmpT = [];
    for idx = 1:1:sizL;
        % Steepness of upper slope [dB/Bark] in accordance with Terhardt
        steep = -24-(230/freqs(idx))+(0.2*LdB(idx_L(idx)));
        steepT = -24-(230/freqs(idx_L(idx)))+(0.2*LdB(idx_L(idx)));
        
        tmpF = [tmpF [freqs(idx);steep]];
        tmpT = [tmpT [freqs(idx_L(idx));steepT]];
        
        if idx == 1     warning('to check formula');    end
        
        if steep < 0
            S2(idx) = steep;
        end
    end

    whichZ	= zeros(2,sizL); % memory allocation
    ExcAmp = zeros(sizL,47); % memory allocation
    Slopes = zeros(sizL,47); % memory allocation
    
    qd		= 1:1:sizL;
    whichZ(1,:)	= floor(2*Barkno(idx_L(qd)+N01)); % in which critical band each level above threshold is
    whichZ(2,:)	= ceil(2*Barkno(idx_L(qd)+N01));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each freq component with a level above threshold. One excitation 
    % pattern for each frequency component:
    for k=1:1:sizL 
        Ltmp = LdB(idx_L(k)); 
        Btmp = Barkno(idx_L(k)+N01);

        for l = 1:1:whichZ(1,k) % lower slope S1, from 1 to lower limit of CBk
            Stemp = (S1*(Btmp-(l*0.5)))+Ltmp;
            if Stemp>MinBf(l)
                Slopes(k,l)=From_dB(Stemp);
            end
        end
        
        for l = whichZ(2,k):1:47 % higher slope S2, from higher limit of CBk to 47
            Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp;
            if Stemp>MinBf(l)
                Slopes(k,l)=From_dB(Stemp);
            end
        end
        % figure; plot((1:47)*.5, To_dB(Slopes(k,:))); xlabel('Critical-band rate [Bark]'); ylabel('Level [dB]'); title(num2str(k))
        disp('')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:1:47 % each critical band k

        etmp    = zeros(1,N);
        etmpExc = zeros(1,N);
        
        for l=1:1:sizL % each level above threshold

            N1tmp = idx_L(l);
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
            etmpExc(N2tmp) = ExcAmp(N1tmp,k);
        end

        if k == 17
            if bDebug
                figure(2); 
                plot(freqs,20*log10(abs(etmp(qb)))) % figure; plot(freqs,20*log10(abs(etmpExc(qb)))), xlabel('Excitation pattern')
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
