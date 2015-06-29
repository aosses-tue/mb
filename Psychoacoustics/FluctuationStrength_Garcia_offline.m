function [dataOut out] = FluctuationStrength_Garcia_offline(insig,fs,N)
% function [dataOut out] = FluctuationStrength_Garcia_offline(insig,fs,N)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
%       Comments on this implementation: cross-correlation seems not to be 
%       working properly
%
%       Comments:
%           fs  - should be an input parameter
%           N   - should be an input parameter
%           Too difficult to change N!
% 
% 2. Stand-alone example:
%       r20141126_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Rodrigo Garcia L./Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 29/06/2015
% Last update on: 29/06/2015 % Update this date manually
% Last use on   : 29/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General fluctuation strength parameters
params  = FluctuationStrength_Garcia_getparams(N);

if nargin < 3
    N   = params.N;
end

% Multiframe separation
overlap = 0.50;
b       = buffer(insig,N,N * overlap,'nodelay');
nFrames = size(b,2);
    
% Blackman window
window = blackman(N);
ampCal = From_dB(100) * 2 / (N * mean(window));
window = ampCal * window';
    
FS = zeros(1,nFrames);
for iFrame = 1:nFrames
    x = b(:,iFrame);

    % Excitation patterns
    tempIn  = window .* x';
    tempIn  = params.a0 .* fft(tempIn);
    Lg      = abs(tempIn(params.qb));
    LdB     = To_dB(Lg);
    whichL  = find(LdB > params.MinExcdB);
    nL      = length(whichL);

    % Steepness of slopes (Terhardt)
    S1 = -27;			
    S2 = zeros(1,nL);
    for w = 1:1:nL;
        steep = -24 - (230 / params.freqs(whichL(w))) ...
            + (0.2 * LdB(whichL(w)));
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

    ei      = zeros(47,N);
    Fei     = zeros(47,N);
    h0      = zeros(1,47);
    hBPi    = zeros(params.Chno,N);
    hBPrms  = zeros(1,params.Chno);
    mdept	= zeros(1,params.Chno);
    Hweight = Test_Hweight;

    for k = 1:1:47
        etmp = zeros(1,N);

        for l = 1:1:nL
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

        ei(k,:)		= N * real(ifft(etmp));
        etmp		= abs(ei(k,:));
        h0(k)		= mean(etmp);
        Fei(k,:)	= fft(etmp - h0(k));
        hBPi(k,:)   = 2 * real(ifft(Fei(k,:) .* Hweight));
        hBPrms(k)   = rmsDik(hBPi(k,:));

        if h0(k) > 0
            mdept(k) =	hBPrms(k) / h0(k);
            if mdept(k) > 1
                mdept(k) = 1;
            end
        else
            mdept(k)=0;
        end
    end

    ki = zeros(1,params.Chno - 2);
    ri = zeros(1,params.Chno);

    % Find cross-correlation coefficients
    for k=1:1:45
        cfac    = cov(hBPi(k,:),hBPi(k + 2,:));
        den     = diag(cfac);
        den     = sqrt(den * den');

        if den(2,1) > 0
            ki(k) = cfac(2,1) / den(2,1);
        else
            ki(k) = 0;
        end
    end

    % Uses test gzi parameter
    gzi = Test_gzi;

    % Calculate specific roughness ri and total roughness R
    ri(1) = (gzi(1) * mdept(1) * ki(1)) ^ 2;
    ri(2) = (gzi(2) * mdept(2) * ki(2)) ^ 2;

    for k = 3:1:45
        ri(k) = (gzi(k) * mdept(k) * ki(k - 2) * ki(k)) ^ 2;
    end

    ri(46) = (gzi(46) * mdept(46) * ki(44)) ^ 2;
    ri(47) = (gzi(47) * mdept(47) * ki(45)) ^ 2;

    FS(iFrame) = params.Cal * sum(ri);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
