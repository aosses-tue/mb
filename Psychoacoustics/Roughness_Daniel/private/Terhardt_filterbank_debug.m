function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank_debug(insig_f,fs,N,qb,MinExcdB,Barkno,N01,zb,Chno,bDebug)
% function [ei, etmp_fd etmpExc_fd] = Terhardt_filterbank_debug(insig_f,Fs,N,qb,MinExcdB,Barkno,N01,zb,Chno,bDebug)
%
% 1. Description:
%       insig_f - input signal in frequency domain.
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       References: Terhardt1979
%       See also: Roughness_offline.m, Terhardt_filterbank.m
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Version       : 1 (based on the same version number as Terhardt_filterbank.m)
% Created on    : 11/02/2016
% Last update on: 11/02/2016 
% Last use on   : 11/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10
    bDebug = 1;
end

if nargin < 9
    Chno = 47;
end

freqs	= (qb+1)*fs/N;
MinBf   = MinExcdB(zb);

Lg		= abs(insig_f);     % It takes only the frequencies of interest
                            % semilogx(freqs,Lg), xlabel('Frequency [Hz]'), ylabel('Magnitude')
Lg      = Lg(qb);
LdB		= To_dB(Lg);
LdB_idx	= find(LdB>MinExcdB); % index of levels above threshold
nL	    = length(LdB_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDebug
    
    hFig = figure(2);
    
    subplot(2,1,1)
    semilogx(freqs,MinExcdB), hold on, grid on
    plot(freqs(LdB_idx),LdB(LdB_idx),'r')
    
    xlim([minmax(freqs(LdB_idx))].*[0.5 2]) % One octave below and above the frequencies above Thres.
    ylim([min(MinExcdB) max(max(LdB),max(MinExcdB))*2])
    xlabel('Frequency [Hz]')
    ylabel('Magnitude')
    
    legend('Hearing threshold','components above threshold')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assessment of slopes (Terhardt) for freq components above thres.
S1 = -27;
S2 = zeros(1,nL);

for idx = 1:1:nL;
    % Steepness of upper slope [dB/Bark] in accordance with Terhardt
    steepT = -24-(230/freqs(LdB_idx(idx)))+(0.2*LdB(LdB_idx(idx)));

    if steepT < 0
        S2(idx) = steepT;
        % if idx == 1     warning('Correction introduced by AO/RG. steep replaced by steepT (see previous control versions of this script');    end
    end
end

whichZ      = zeros(2,nL)   ;   % memory allocation
ExcAmp      = zeros(nL,Chno);   % memory allocation
Slopes      = zeros(nL,Chno);   % memory allocation
etmp_fd     = zeros(Chno,N) ;
etmpExc_fd  = zeros(Chno,N) ;

qd          = 1:1:nL;
% In which critical band number are the levels above threshold located:
whichZ(1,:)	= floor(2*Barkno(LdB_idx(qd)+N01)); % idxs up to which S1 is going to be used
whichZ(2,:)	= ceil(2*Barkno(LdB_idx(qd)+N01)); % idxs from which S2 is going to be used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each freq component with a level above threshold. One excitation 
% pattern for each frequency component. Therefore sizL patterns will be obtained.
% The masking patterns are going to be in the critical bands between
% whichZ(1,:) and whichZ(2,:)
for k=1:1:nL 
    Ltmp = LdB(LdB_idx(k));         % Level [dB]
    Btmp = Barkno(LdB_idx(k)+N01);  % Frequency [Bark]

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

if bDebug
    
    figure(2);
    subplot(2,1,2)
    
    idxq = LdB_idx;
    L = length(idxq);
    
    idx2plot = 14:22;
    BandNumber = repmat( ((idx2plot)'),1,L);
    
    lvl_idx = 1:length(freqs(LdB_idx));
    
    level_sm = repmat(lvl_idx,length(idx2plot),1);
    mesh(level_sm',BandNumber',To_dB( Slopes(:,idx2plot) )), hold on
    
    xlim(minmax(lvl_idx))
    
    xlabel(sprintf('Level_{idx}\nnumber of levels above threshold = %.0f',max(lvl_idx)))
    ylabel('CB number')
    title('Determined masking patterns')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:Chno % each critical band k

    etmp    = zeros(1,N);
    etmpExc = zeros(1,N);

    for l=1:1:nL % each level above threshold

        N1tmp = LdB_idx(l);
        N2tmp = N1tmp + N01;
        
        if (whichZ(1,l) == k)
            
            ExcAmp(N1tmp, k) = 1; % Increasing slopes, increasing Lg as k increases:
                                  % So: ExpAmp increases with values from 0 to 1
            if bDebug; idxExcAmp(l,k)=1; end
            
        elseif (whichZ(2,l) == k) % If lower limit is equal to k (Second enters here)
            
            ExcAmp(N1tmp, k) = 1;
            if bDebug; idxExcAmp(l,k)=2; end
            
        elseif (whichZ(2,l) > k) % If upper limit is equal to k (Third enters here)
            
            ExcAmp(N1tmp,k) = Slopes(l,k+1)/Lg(N1tmp);
            if bDebug; idxExcAmp(l,k)=3; end
            
        else
            
            ExcAmp(N1tmp,k) = Slopes(l,k-1)/Lg(N1tmp); % Decreasing slopes, decreasing Lg as k decreases
                                                       % So: ExpAmp decreases with values from 1 to 0
            if bDebug; idxExcAmp(l,k)=4; end
            
        end
        etmp(N2tmp)     = ExcAmp(N1tmp,k)*insig_f(N2tmp); 
        etmpExc(N2tmp)  = ExcAmp(N1tmp,k);
    end

    if bDebug; disp(sprintf('CB number k = %.0f, ExcAmp=%.2f, entered ''if'' %.0f',k,ExcAmp(l, k),idxExcAmp(l,k))); end
        
    % figure; ktmp = 13; plot(Slopes(:,ktmp)./transpose(Lg(whichL)));
    % figure; plot(freqs, abs(etmp(qb)) ), hold on; plot(freqs(whichL),ExcAmp(:,ktmp)/max(ExcAmp(:,ktmp))*max(abs(etmp)),'r'), xlim(minmax(freqs(whichL)).*[0.8 1.2]), pause(2), close
    
    if k == 17
        if bDebug
            hFig(end+1) = figure(3); 
            plot(freqs,20*log10(abs(etmp(qb))),'LineWidth',2), hold on
            plot(freqs,20*log10(abs(insig_f(qb))),'r--')
            hold on; grid on
            xlabel('Frequency [Hz]')
            ylabel('Level [dB]')
            title('Spectral components: Terhardt''s model')
            legend('Excitation, etmp (to be used for calculation)','M + all freq components')
            xlim([900-200 1100+200])
        end
    end
    
    etmp_fd(k,:)    = etmp;
    etmpExc_fd(k,:) = etmpExc;
    ei(k,:)     = N*real(ifft( etmp_fd(k,:) )); % N factor as in the Dik's code
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
