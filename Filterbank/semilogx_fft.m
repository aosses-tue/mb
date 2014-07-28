function [H_dB, freqs] = semilogx_fft(B, A, fs, fftsize, bPlot, bPhase)
% function [H_dB, freqs] = semilogx_fft(B, A, fs, fftsize, bPlot, bPhase)
% 
% Plots FFT of an IIR filter defined by B (num coeffs), A (den coeffs), 
% discarding the redundant values, i.e., it plots first K bins where K = fftsize/2+1
%
% % Example:
%       [x, Fs] = wavread('/home/alejandro/Documenten/Meas/Meas/Music/Loudness_balanced_stimuli/PR_Stimuli_SJ/UW_LB_F0m_131_Hz.wav');
%       plotfft_from_coeff(x,1,4096,Fs,0)
%
% % Standalone example:
%       [x, Fs] = wavread('/home/alejandro/Documenten/Meas/Meas/Music/Loudness_balanced_stimuli/PR_Stimuli_SJ/UW_LB_F0m_131_Hz.wav');
%       [H_dB, f] = semilogx_fft(x,1,Fs);
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return; % Just testing from Check_dependencies_ExpORL
end

if nargin   < 6
    bPhase  = 0;
end

if nargin   < 5
    bPlot  = 0;
end

if nargin   < 4
    fftsize = 4096*2;
end

if nargin   < 3
    fs      = 44100;
end

if nargin   < 2
    A       = 1;
end

freqs   = [1:fftsize/2]'/fftsize*fs;
[size_row  size_col ] = size(B);
[size_rowA size_colA] = size(A);

if size_colA ~= size_col
    A = repmat(A,1,size_col);
end


for i = 1:size_col 
    h(:,i)  = freqz(B(:,i), A(:,i), fftsize/2);
    phase   = unwrap(angle(h(i)));
end

if length(phase == 1) %linear phase
    phase = repmat(phase,length(h),1);
end

delay   = diff(phase)*180/pi./diff(freqs);

if bPhase == 1 && bPlot == 1
    figure
    a(1)    = subplot(2,1,1);
end

H_dB = 20*log10(abs(h));

if bPlot == 1
    semilogx(freqs, H_dB);
    ylabel('Magnitude [dB]');
    grid on
end

if bPhase == 1
    a(2)    = subplot(2,1,2);
    plotyy(freqs, phase, freqs, [0; delay]);
    grid on
    linkaxes(a,'x');
end
