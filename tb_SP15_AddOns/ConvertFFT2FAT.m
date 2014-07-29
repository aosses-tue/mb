function [freq_warped, H_dB, freq] = ConvertFFT2FAT(filename)
% function [freq_warped, H_dB, freq] = ConvertFFT2FAT(filename)
%
% 1. Description:
%       freq is mapped into freq_warped. freq_warped is linearly interpolated
%       between fi and fs of each FAT-channel
% 
% Programmed by Alejandro Osses
% Created in: 2013
% Last updated on: 21/05/2014
% Last used on   : 21/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    filename = '/home/alejandro/Documenten/Meas/Meas/Music/Loudness_balanced_stimuli/PR_Stimuli_SJ/UW_LB_F0m_131_Hz.wav';
end

[x, Fs] = wavread(filename);
N = 2^16;
[H_dB, freq] = semilogx_fft(x,1,Fs,N);

p.CIC.FAT   = 22; % Number of channels
p.CFG.Fs    = 15659.375; % fs xPC Target

FAT         = nsb_ReadFAT(p.CIC.FAT).*p.CFG.Fs./16000;

freq_warped = [];
disc_begin  = find(freq<FAT(1));
disc_end    = find(freq>FAT(23));

for i = 1:length(FAT)-1
    idx = find(freq>=FAT(i) & freq<=FAT(i+1));
    freq_warped = [freq_warped; ( freq(idx)-freq(min(idx)) )./(freq(max(idx))-freq(min(idx)))+i];
end

H_dB(disc_end)   = [];
H_dB(disc_begin) = [];

H_dB    = H_dB - max(H_dB);
idx     = find(H_dB<-40);
H_dB(idx) = -40;

H_dB    = ( H_dB + 40 )/40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end