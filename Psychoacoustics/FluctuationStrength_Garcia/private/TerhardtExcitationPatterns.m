function ei = TerhardtExcitationPatterns(insig,fs,dBFS)
% function ei = TerhardtExcitationPatterns(insig,fs,dBFS)
% 
% 1. Description:
%       Calculates Terhardt excitation patterns.
%       Proper scaling:
%           If the input signal is a AM sine tone (fc = 1 kHz, fmod = 4 Hz, 
%           full modulation) at a level of 60 dB the three spectral components
%           will be at 1000 Hz (58.2 dB) and the side components at 996 and 
%           1004 Hz respectively (52.3 dB each). Then the total level can 
%           be re-calculated as sum_db([58.2 52.3 52.3]). You can use an AM tone
%           to confirm that the scaling is correct.
% 
% Inputs:
%   insig: The input signal.
% 
% Outputs:
%   ei: Excitation patterns.
% 
% Author: Alejandro Osses Vecchi
% version: 2 on 5/02/2016 (not updated in general script +3 in line 24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dBFS = 100; % AMT toolbox convention
end
corr = dBFS + 3;

% General parameters
params = il_calculate_params(insig,fs);

% Transforms input signal to frequency domain
insig = From_dB(corr)*fft(insig)/params.N; % 3 dB added to adjust the SPL values to be put into slope equations

% Use only samples that fall into the audible range
Lg  = abs(insig(params.qb));
LdB = To_dB(Lg);

% Use only components that are above the hearing threshold
whichL = find(LdB > params.MinExcdB);
nL     = length(whichL);

% Steepness of slopes
S1 = -27;			
S2 = zeros(1,nL);
for w = 1:nL;
    steep = -24 - (230 / params.freqs(whichL(w))) ...
        + (0.2 * LdB(whichL(w)));
    if steep < 0
        S2(w) = steep;
    end
end

whichZ      = zeros(2,nL);
whichZ(1,:)	= floor(2 * params.Barkno(whichL + params.N01));
whichZ(2,:)	= ceil(2 * params.Barkno(whichL + params.N01));

% Calculate slopes from steep values
Slopes = zeros(nL,params.Chno);
for l = 1:nL
    Ltmp = LdB(whichL(l));
    Btmp = params.Barkno(whichL(l) + params.N01);

    for k = 1:whichZ(1,l)
        Stemp =	(S1 * (Btmp - (k * 0.5))) + Ltmp;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = From_dB(Stemp);
        end
    end

    for k = whichZ(2,l):params.Chno
        Stemp = (S2(l) * ((k * 0.5) - Btmp)) + Ltmp;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = From_dB(Stemp);
        end
    end 
end

% Excitation patterns
ExcAmp  = zeros(nL,params.Chno);
ei      = zeros(params.Chno,params.N);
for k = 1:params.Chno
    etmp = zeros(1,params.N);
    for l = 1:nL
        N1tmp = whichL(l);
        N2tmp = N1tmp + params.N01;

        if whichZ(1,l) == k
            ExcAmp(N1tmp,k)	= 1;
        elseif whichZ(2,l) == k
            ExcAmp(N1tmp,k)	= 1;
        elseif whichZ(2,l) > k
            ExcAmp(N1tmp,k) = Slopes(l,k+1)/Lg(N1tmp);
        else % whichZ(1,l) < k
            ExcAmp(N1tmp,k) = Slopes(l,k-1)/Lg(N1tmp);
        end

        etmp(N2tmp) = ExcAmp(N1tmp,k)*insig(N2tmp);
    end

    ei(k,:) = 2*params.N*real(ifft(etmp));
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = il_calculate_params(x,fs)
    params      = struct;
    params.N    = length(x);
    params.Chno = 47;

    % Defines audible range indexes and frequencies
    df           = fs/params.N;
    N0           = round(20/df)+1;
    Ntop         = round(20e3/df)+1;
    params.N01   = N0-1;
    params.qb    = N0:Ntop;
    params.freqs = (params.qb-1)*df;

    [params.Barkno Bark_raw] = Get_Bark(params.N,params.qb,params.freqs);
    % Loudness threshold related parameters
    params.MinExcdB = calculate_MinExcdB(params.N01,params.qb,params.Barkno);
    params.MinBf    = calculate_MinBf(params.N01,df,Bark_raw,params.MinExcdB);
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinExcdB = calculate_MinExcdB(N01,qb,Barkno)
    HTres = [
        0		130
        0.01    70
        0.17    60
        0.8     30
        1       25
        1.5     20
        2		15
        3.3     10
        4		8.1
        5		6.3
        6		5
        8		3.5
        10		2.5
        12		1.7
        13.3	0
        15		-2.5
        16		-4
        17		-3.7
        18		-1.5
        19		1.4
        20		3.8
        21		5
        22		7.5
        23      15
        24      48
        24.5 	60
        25		130
    ];

    MinExcdB            = zeros(1,length(qb));
    MinExcdB(qb-N01)    = interp1(HTres(:,1),HTres(:,2),Barkno(qb));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinBf = calculate_MinBf(N01,df,Bark,MinExcdB)
    Cf = round(Bark(2:25,2)'/df)-N01+1;
    Bf = round(Bark(1:25,3)'/df)-N01+1;  

    zb      = sort([Bf Cf]);
    MinBf   = MinExcdB(zb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
