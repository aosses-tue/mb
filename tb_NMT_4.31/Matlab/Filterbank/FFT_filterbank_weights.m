function weights = FFT_filterbank_weights(channel_spec, fft_length, audio_sample_rate)

% FFT_filterbank_weights: Calculate FFT filterbank weights
%
% weights = FFT_filterbank_weights(channel_bounds, fft_length, audio_sample_rate)
% weights = FFT_filterbank_weights(channel_edges,  fft_length, audio_sample_rate)
%
% fft_length:        FFT length.
% audio_sample_rate: Audio sample rate (in Hertz).
% channel_bounds:    Frequency boundaries of filter bands; matrix size = (num_channels, 2)
%                    Each row is a pair of frequencies: the upper and lower boundary.
% channel_edges:     Frequency boundaries of filter bands; matrix size = (num_channels+1, 1)
% weights:           Weights matrix used to combine FFT bins into filter bands.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_bounds = FFT_bin_bounds(fft_length, audio_sample_rate);
[num_bins,     num_cols] = size(bin_bounds);
if num_cols ~= 2
	error('bin_bounds must have two columns');
end
bin_width = audio_sample_rate / fft_length;

[num_rows, num_cols] = size(channel_spec);
switch num_cols
case 1
	% Edge frequencies specified:
	channel_bounds = [channel_spec(1:end-1), channel_spec(2:end)];
	num_channels = num_rows - 1;
case 2
	% Upper & lower boundary frequencies specified:
	channel_bounds = channel_spec;
	num_channels = num_rows;
otherwise
	error('channel_bounds must have two columns');
end

weights = zeros(num_channels, num_bins);

for chan = 1:num_channels
	for bin = 1:num_bins
		% Lower bound is the maximum of channel lower bound and bin lower bound:
		f_lower = max(bin_bounds(bin, 1), channel_bounds(chan, 1));
		% Upper bound is the minimum of channel upper bound and bin upper bound:
		f_upper = min(bin_bounds(bin, 2), channel_bounds(chan, 2));
		
		% Weight is the proportion of the bin that is spanned by the channel:
        f_width = f_upper - f_lower;
		if (f_width > 0)
			weights(chan, bin) = f_width / bin_width;
		end
	end
end
