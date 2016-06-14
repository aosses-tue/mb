function outsig = PeripheralHearingSystem(insig,fs,N)
% function y = PeripheralHearingSystem(insig,fs,N)
% 
% Applies the effect of transmission from free field to the cochlea to a
% given signal.
% 
% Inputs:
% insig: The signal to process.
% fs: Sampling frequency,
% N: number of samples.
% 
% Author: Rodrigo García

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0 = il_calculate_a0(fs,N);
outsig  = 2*real(ifft(a0.*fft(insig)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a0 = il_calculate_a0(fs,N)
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
