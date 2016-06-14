function outsig = PeripheralHearingSystem_t(insig,fs,bIdle)
% function y = PeripheralHearingSystem_t(insig,fs,bIdle)
% 
% Applies the effect of transmission from free field to the cochlea to a
% given signal. Time domain version.
% 
% Inputs:
% insig: The signal to process. insig has to be a row vector.
% fs: Sampling frequency,
% 
% Author: Alejandro Osses V.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 
    bIdle = 0;
end

K = 2^12;

if bIdle == 0; % normal processing including resonance of the ear canal
    B = il_calculate_a0(fs,K);
end

if bIdle == 1;
    B = il_calculate_a0_idle(fs,K);
end

outsig = filter(B,1,[insig zeros(1,K/2)]);
outsig = outsig(K/2+1:end);

if nargout == 0
    if bIdle == 0
        il_calculate_a0(fs,K);
    else
        il_calculate_a0_idle(fs,K);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_calculate_a0(fs,N)
% Compensation of the transmission factor from Free-field, taken from
% Fastl2007, Fig. 8.18, page 226

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = Get_Bark(N,qb,freqs);

a0tab = [
    0       0
    10      0
    12      1.15
    13      2.31
    14      3.85
    15      5.62
    16      6.92
    16.5    7.38
    17      6.92
    18      4.23
    18.5    2.31
    19      0
    20      -1.43
    21		-2.59
    21.5	-3.57
    22		-5.19
    22.5	-7.41
    23		-11.3
    23.5	-20
    24		-40
    25		-130
    26		-999
];

a0            = zeros(1,N);
a0(qb)        = From_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = il_create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    il_create_a0_FIR(freqs,a0(qb),N,fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_calculate_a0_idle(fs,N)
% No resonance of the ear canal accounted for.

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = Get_Bark(N,qb,freqs);

a0tab = [
    0       0
    10      0
    19      0
    20      -1.43
    21		-2.59
    21.5	-3.57
    22		-5.19
    22.5	-7.41
    23		-11.3
    23.5	-20
    24		-40
    25		-130
    26		-999
];

a0            = zeros(1,N);
a0(qb)        = From_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = il_create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    il_create_a0_FIR(freqs,a0(qb),N,fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_create_a0_FIR(f,a0,N,fs)

f = [0 f fs/2];
a0 = [a0(1) a0 a0(end)];

B = fir2(N,f/(fs/2),a0);

if nargout == 0
    [H1,Fn]=freqz(B,1,N/2);
    
    figure;
    plot(fs/2*Fn/pi, 20*log10(abs([H1])));
    xlabel('Frequency [Hz]')
    legend([num2str(N) ' taps']);
    title('FIR filter to be used as approximation to isolation curve')
    xlim([0 fs/2])
end
