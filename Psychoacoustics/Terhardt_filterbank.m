function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank(TempIn,Fs,N,qb,MinExcdB,Barkno,N01,zb)
% function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank(TempIn,Fs,N,qb,MinExcdB,Barkno,N01,zb)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 26/05/2015
% Last update on: 26/05/2015 % Update this date manually
% Last use on   : 25/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freqs	= (qb+1)*Fs/N;
MinBf   = MinExcdB(zb);

Lg		= abs(TempIn(qb));
LdB		= To_dB(Lg);
idx_L	= find(LdB>MinExcdB); % index of levels above threshold
sizL	= length(idx_L);

% steepness of slopes (Terhardt)
S1 = -27;
S2 = zeros(1,sizL);

tmpF = [];
tmpT = [];
for idx = 1:1:sizL;
    % Steepness of upper slope [dB/Bark] in accordance with Terhardt
    steep   = -24-(230/freqs(idx))       +(0.2*LdB(idx_L(idx))); % original from Dik's code
    steepT  = -24-(230/freqs(idx_L(idx)))+(0.2*LdB(idx_L(idx)));

    tmpF = [tmpF [freqs(idx);steep]];
    tmpT = [tmpT [freqs(idx_L(idx));steepT]];

    if steepT < 0
        S2(idx) = steepT;
        if idx == 1     warning('Correction introduced by AO/RG. steep replaced by steepT');    end
    end
end

whichZ      = zeros(2,sizL);    % memory allocation
ExcAmp      = zeros(sizL,47);   % memory allocation
Slopes      = zeros(sizL,47);   % memory allocation
etmp_fd     = zeros(47,N);
etmpExc_fd  = zeros(47,N);

qd          = 1:1:sizL;
whichZ(1,:)	= floor(2*Barkno(idx_L(qd)+N01)); % in which critical band each level above threshold is
whichZ(2,:)	= ceil(2*Barkno(idx_L(qd)+N01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if k == 17
        disp('')
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        etmp(N2tmp)     = ExcAmp(N1tmp,k)*TempIn(N2tmp); 
        etmpExc(N2tmp)  = ExcAmp(N1tmp,k);
    end

    etmp_fd(k,:)    = etmp;
    etmpExc_fd(k,:) = etmpExc;
    ei(k,:)     = N*real(ifft( etmp_fd(k,:) ));
    
    %% end Terhardt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
