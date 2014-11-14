function bin_bounds = FFT_bin_bounds(fft_length, audio_sample_rate)

% FFT_bin_bounds: Calculate FFT bin frequency boundaries
%
% bin_bounds = FFT_bin_bounds(fft_length, audio_sample_rate)
%
% fft_length:        FFT length.
% audio_sample_rate: Audio sample rate (in Hertz).
% bin_bounds:        Frequency boundaries of FFT bins; matrix size = (num_bins, 2)
%
% Each row of the bin_bounds matrix is a pair of frequencies: upper and lower boundary.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_spacing = audio_sample_rate/fft_length;
K = fft_length/2;
k = (0:K)';
kk = [k - 0.5, k + 0.5];
bin_bounds = kk * freq_spacing;
