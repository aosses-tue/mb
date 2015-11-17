function [dataOut out] = FluctuationStrength_Garcia_offline(insig, fs, N, optsDebug)
% function [dataOut out] = FluctuationStrength_Garcia_offline(insig, fs, N, optsDebug)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
% 
%       Comments on this implementation: cross-correlation seems not to be 
%       working properly.
%
%       Adapted from Model/Helper/FluctuationStrength.m
% 
%       Comments:
%           fs  - should be an input parameter
%           N   - should be an input parameter
%           Too difficult to change N!
%           gzi - where was it taken from?
% 
% 2. Stand-alone example:
%       r20141126_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Rodrigo Garcia L./Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 29/06/2015
% Last update on: 16/11/2015 
% Last use on   : 16/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    optsDebug = [];
end

optsDebug   = ef(optsDebug,'all',0);
optsDebug   = ef(optsDebug,'ki',0);

bDebug      = optsDebug.all;
debug_ki    = optsDebug.ki;

if optsDebug.all == 1 | optsDebug.ki == 1
    disp('Debug mode active. Only the first frame will be analysed');
    hFig = [];
end

% General fluctuation strength parameters
params  = FluctuationStrength_Garcia_getparams(N);

if nargin < 3
    N   = params.N;
end

Chno = params.Chno;

% Multiframe separation
overlap = 0.50;
b       = buffer(insig,N,N * overlap,'nodelay');
nFrames = size(b,2);
    
% Blackman window
window  = blackman(N);
ampCal  = From_dB(100) * 2 / (N * mean(window));
window  = ampCal * window';
    
FS      = zeros(1,nFrames);

for iFrame = 1:nFrames
    
    h0      = zeros(   1,Chno);
    hBPi    = zeros(Chno,   N);
    mdept	= zeros(   1,Chno);
    
    x = b(:,iFrame);

    SPL(iFrame) = rmsdb(x)+100;
    tempIn  = window .* x';
    
    %% Peripheral stage:
    
    % Transmission factor a0:
    tempIn  = params.a0 .* fft(tempIn);
    
    % Excitation patterns:
    Lg      = abs(tempIn(params.qb));
    LdB     = To_dB(Lg);
    whichL  = find(LdB > params.MinExcdB);
    nL      = length(whichL);

    % Steepness of slopes (Terhardt)
    S1 = -27;			
    S2 = zeros(1,nL);
    for w = 1:1:nL;
        steep = -24 - (230 / params.freqs(whichL(w))) + (0.2 * LdB(whichL(w)));
        if steep < 0
            S2(w) = steep;
        end
    end

    whichZ      = zeros(2,nL);
    qd          = 1:1:nL;
    whichZ(1,:)	= floor(2 * params.Barkno(whichL(qd) + params.N01));
    whichZ(2,:)	= ceil(2 * params.Barkno(whichL(qd) + params.N01));

    ExcAmp = zeros(nL,47);
    Slopes = zeros(nL,47);

    %%%
    for k = 1:1:nL    
        Ltmp = LdB(whichL(k));
        Btmp = params.Barkno(whichL(k) + params.N01);

        for l = 1:1:whichZ(1,k)
            Stemp =	(S1 * (Btmp - (l * 0.5))) + Ltmp;
            if Stemp > params.MinBf(l)
                Slopes(k,l) = From_dB(Stemp);
            end
        end

        for l = whichZ(2,k):1:47
            Stemp = (S2(k) * ((l * 0.5) - Btmp)) + Ltmp;
            if Stemp > params.MinBf(l)
                Slopes(k,l) = From_dB(Stemp);
            end
        end 
    end
    %%% 
    
    for k = 1:Chno
        etmp = zeros(1,N);

        for l = 1:nL
            N1tmp = whichL(l);
            N2tmp = N1tmp + params.N01;

            if whichZ(1,l) == k
                ExcAmp(N1tmp,k)	= 1;
            elseif whichZ(2,l) == k
                ExcAmp(N1tmp,k)	= 1;
            elseif whichZ(2,l) > k
                ExcAmp(N1tmp,k) = Slopes(l,k + 1) / Lg(N1tmp);
            else
                ExcAmp(N1tmp,k) = Slopes(l,k - 1) / Lg(N1tmp);
            end

            etmp(N2tmp) = ExcAmp(N1tmp,k) * tempIn(N2tmp);
        end

        if bDebug
            etmp_fd(k,:)= etmp;
        end
        etmp_td(k,:)= N * real(ifft(etmp)); %etmp_td - output of the critical-band filter bank
    end
    
    if bDebug == 0
        clear ExcAmp Slopes; % it liberates some memory
    end
    %%% 
    
    Hweight = il_create_Hweight(N,fs);

    for k = 1:Chno
        
        h0(k)		= mean(abs(etmp_td(k,:)));
        Fei         = fft(abs(etmp_td(k,:)) - h0(k));
        hBPi(k,:)   = 2 * real(ifft(Fei.* Hweight'));
        if bDebug
            hBPino(k,:) = 2 * real(ifft(Fei       ));
        end
        
    end
    clear Hweight etmp;
    
    hBPrms  = rms(transpose(hBPi));
    if bDebug
        hBPrmsno= rms(transpose(hBPino));
    end
    
    idx = find(h0 > 0);
    mdept(idx) = hBPrms(idx) ./ h0(idx);
    if bDebug
        mdeptno(idx) = hBPrmsno(idx) ./ h0(idx);
    end
    
    idx = find(mdept > 1);
    mdept(idx) = 1;
    
    if bDebug
    
        figM = 3;
        figN = 2;
        plot_x = 1:size(etmp_td,2);
        plot_f = ( 1:size(etmp_td,2) )/N*fs;

        %%%
        for k = 1:1:Chno
        
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
                grid on
                ha(end+1) = gca;
                linkaxes(ha,'xy');
                ylim([18 62])
                xlim([3000/N*fs 9000/N*fs]) % plotted from bin 3000 to 9000
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
        %%%
               
    end
    
    ki = zeros(1,Chno - 2);
    fi = zeros(1,Chno);

    % Find cross-correlation coefficients
    for k=1:1:Chno-2
        cfac    = cov(hBPi(k,:),hBPi(k + 2,:));
        den     = diag(cfac);
        den     = Round( sqrt(den * den'), 10); % rounded to 10 decimals

        if den(2,1) > 0
            ki(k) = cfac(2,1) / den(2,1);
        elseif den(2,1) < 0
            ki(k) = 0;
        else
            ki(k) = 0;
        end
    end

    if debug_ki & iFrame == 1
    
        ki_fac( 1: 2) = ki(1:2).^2;
        ki_fac( 3:45) = ki(3:45).*ki(1:43);
        ki_fac(46:47) = ki(44:45).^2;
        
        zi = params.zi;
        
        figure;
        plot(zi,ki_fac,'LineWidth',2); hold on
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
    
    % Uses test gzi parameter
    % gzi = Test_gzi(N);
    gzi = il_create_gzi(Chno);
    
    % Calculate specific roughness ri and total roughness R
    fi(iFrame,1) = (gzi(1) * mdept(1) * ki(1)) ^ 2;
    fi(iFrame,2) = (gzi(2) * mdept(2) * ki(2)) ^ 2;

    for k = 3:1:45
        fi(iFrame,k) = (gzi(k) * mdept(k) * ki(k - 2) * ki(k)) ^ 2;
    end

    fi(iFrame,46) = (gzi(46) * mdept(46) * ki(44)) ^ 2;
    fi(iFrame,47) = (gzi(47) * mdept(47) * ki(45)) ^ 2;

    FS(iFrame) = params.Cal * sum(fi(iFrame,:));
end

dataOut{1} = FS;
dataOut{2} = fi;
dataOut{3} = SPL;

nParam      = 1;
out.Data1   = transpose(FS);
output.name{nParam} = 'Fluctuation strength';
output.param{nParam} = strrep( lower( output.name{nParam} ),' ','-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

function Hweight = il_create_Hweight(N,fs)

params = struct;
% params.sf1 = 0.5;
% params.pf1 = 2; % 2
% params.pf2 = 8; % 16
% params.sf2 = 32; % 32
params.sf1 = 0.5;
params.pf1 = 2; % 2
params.pf2 = 8; % 16
params.sf2 = 64; % 32

d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',params.sf1,params.pf1,params.pf2,params.sf2,80,3,80,fs);
Hd = design(d,'butter');
% fvtool(Hd);
% measure(Hd)

x = [1; zeros(N-1,1)];
y = filter(Hd,x);
freq = (0:(2*pi)/length(x):pi)/pi*fs/2;
xdft = fft(x);
ydft = fft(y);

Hweight = abs(ydft);
HweightdB = To_dB(Hweight); 

warning('not exactly as implemented by Rodrigo (MATLAB 2015 needed');

% figure;
% plot(freq,20*log10(abs(ydft(1:length(x)/2+1))),'r','linewidth',2);
% xlim([0 50])
% disp('')

% legend('Original Signal','Bandpass Signal');

% Hweight = designfilt(...
%         'bandpassiir', ...
%         'StopbandFrequency1', params.sf1, ...
%         'PassbandFrequency1', params.pf1, ...
%         'PassbandFrequency2', params.pf2, ...
%         'StopbandFrequency2', params.sf2, ...
%         'StopbandAttenuation1', 100, ...
%         'PassbandRipple', 3, ...
%         'StopbandAttenuation2', 100, ...
%         'SampleRate', fs);

function gzi = il_create_gzi(Chno);

    Chstep = 0.5;
    
    g0 = [
        0       0
        125     1
        250     1
        500     1
        1000    1
        1500    1
        2000    1
        3000    1
        4000    1
        6000    1
        8000    1
        16000   0
    ];
    
    gzi = interp1(freqtoaud(g0(:,1),'bark'),g0(:,2),(1:Chno)*Chstep);
    gzi(isnan(gzi)) = 0;
    
    gzi = ones(1,Chno);