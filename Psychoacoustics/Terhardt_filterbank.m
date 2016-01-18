function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank(insig_f,fs,N,qb,MinExcdB,Barkno,N01,zb,Chno)
% function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank(insig_f,Fs,N,qb,MinExcdB,Barkno,N01,zb,Chno)
%
% 1. Description:
%       insig_f - input signal in frequency domain.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       References: Terhardt1979
%       See also: Roughness_offline.m
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 26/05/2015
% Last update on: 26/05/2015 
% Last use on   : 06/01/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    Chno = 47;
end

freqs	= (qb+1)*fs/N;
MinBf   = MinExcdB(zb);

Lg		= abs(insig_f(qb));
LdB		= To_dB(Lg);
LdB_idx	= find(LdB>MinExcdB); % index of levels above threshold
nL	    = length(LdB_idx);

% steepness of slopes (Terhardt)
S1 = -27;
S2 = zeros(1,nL);

for idx = 1:1:nL;
    % Steepness of upper slope [dB/Bark] in accordance with Terhardt
    steepT  = -24-(230/freqs(LdB_idx(idx)))+(0.2*LdB(LdB_idx(idx)));

    if steepT < 0
        S2(idx) = steepT;
        if idx == 1     warning('Correction introduced by AO/RG. steep replaced by steepT (see previous control versions of this script');    end
    end
end

whichZ      = zeros(2,nL)   ;    % memory allocation
ExcAmp      = zeros(nL,Chno);   % memory allocation
Slopes      = zeros(nL,Chno);   % memory allocation
etmp_fd     = zeros(Chno,N) ;
etmpExc_fd  = zeros(Chno,N) ;

qd          = 1:1:nL;
whichZ(1,:)	= floor(2*Barkno(LdB_idx(qd)+N01)); % in which critical band each level above threshold is
whichZ(2,:)	= ceil(2*Barkno(LdB_idx(qd)+N01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each freq component with a level above threshold. One excitation 
% pattern for each frequency component:
for k=1:1:nL 
    Ltmp = LdB(LdB_idx(k)); 
    Btmp = Barkno(LdB_idx(k)+N01);

    % Stemp - Masker level L'_{\nu\mu}
    % Anumu - A_{\mu}, Terhardt1979, Eq 6
    for l = 1:1:whichZ(1,k) % lower slope S1, from 1 to lower limit of CBk
        Stemp = (S1*(Btmp-(l*0.5)))+Ltmp; % Terhardt1979, Eq. (4.b)
        if Stemp>MinBf(l)
            Slopes(k,l)=From_dB(Stemp);
        end
    end

    % Stemp - Masker level L''_{\nu\mu}
    for l = whichZ(2,k):1:Chno % higher slope S2, from higher limit of CBk to 47
        Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp; % Terhardt1979, Eq. (4.a)
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
for k=1:1:Chno % each critical band k

    etmp    = zeros(1,N);
    etmpExc = zeros(1,N);

    for l=1:1:nL % each level above threshold

        N1tmp = LdB_idx(l);
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
        etmp(N2tmp)     = ExcAmp(N1tmp,k)*insig_f(N2tmp); 
        etmpExc(N2tmp)  = ExcAmp(N1tmp,k);
    end

    etmp_fd(k,:)    = etmp;
    etmpExc_fd(k,:) = etmpExc;
    ei(k,:)     = N*real(ifft( etmp_fd(k,:) )); % N factor as in the Dik's code
    
    %% end Terhardt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
