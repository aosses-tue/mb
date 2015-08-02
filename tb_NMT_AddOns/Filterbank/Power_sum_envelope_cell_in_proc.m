function out = Power_sum_envelope_cell_in_proc(p, cell_in)
% function out = Power_sum_envelope_cell_in_proc(p, cell_in)
%
% Power-sum envelopes for FFT filterbank.
%
% Inputs:
% p:                Parameter struct containing the fields:
% p.block_length:   The FFT length.
% p.window:         The FFT window function.
% p.num_bins:       The number of FFT bins to retain.
% p.weights:        The weights used in summing the bins.
% p.num_bands:      The number of bands in the filterbank.
% u:                Complex filterbank output FTM
%
% Outputs:
% out{1}:           Real envelope FTM.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from:
%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	out = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fundamental parameters:
    p = Ensure_field(p, 'equalise',  1);
	p = Ensure_field(p, 'num_bands', 22);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Derived parameters:

    if p.block_length == 128;    
        p.band_bins = FFT_band_bins(p.num_bands)';
    elseif p.block_length == 512
        p.band_bins = [2 2 2 3 3 3 4 5 5 6 7 8 10 11 13 15 17 20 23 27 31 35]';
    else
        error('Only block widths of 128 or 512 are allowed.');
    end
	% Weights matrix for combining FFT bins into bands:
	% (incorporates frequency response equalisation)
    
    if ~isfield(p,'weights')
        freq_response  = freqz(p.window/2, 1, p.block_length);
        power_response = freq_response .* conj(freq_response);

        % Check alternative calculation using fft

        P1 = power_response(1);
        P2 = 2 * power_response(2);
        P3 = power_response(1) + 2 * power_response(3);

        p.power_gains = zeros(p.num_bands, 1);
        for band = 1:p.num_bands;
            width = p.band_bins(band);
            if (width == 1)
                p.power_gains(band) = P1;
            elseif (width == 2)
                p.power_gains(band) = P2;
            else
                p.power_gains(band) = P3;
            end;
        end;

        p.weights = zeros(p.num_bands, p.num_bins);
        if p.block_length == 128
            bin = 3;	% We always ignore bins 0 (DC) & 1.
        elseif p.block_length == 512
            bin = 5;
        else
            error('Only block widths of 128 or 512 are allowed.');
        end
        for band = 1:p.num_bands;
            width = p.band_bins(band);
            p.weights(band, bin:(bin + width - 1)) = 1 / p.power_gains(band);
            bin = bin + width;
        end;

        if p.block_length == 128
            cum_num_bins = [1.5; 1.5 + cumsum(p.band_bins)];
        elseif p.block_length == 512
            cum_num_bins = [3.5; 3.5 + cumsum(p.band_bins)];
        else
            error('Only block widths of 128 or 512 are allowed.');
        end

        p.crossover_freqs = cum_num_bins * p.bin_freq;
        p.band_widths = diff(p.crossover_freqs);
        p.char_freqs = p.crossover_freqs(1:p.num_bands) + p.band_widths/2;
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	out = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Essential step for Power_sum_envelope procedure
    u = cell_in{1};
	v = u .* conj(u);						% Power (magnitude squared) of each bin.
	u = p.weights * v;						% Weighted sum of bin powers.
	v = sqrt(u);							% Magnitude.
    out = {v, cell_in{2}};
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;
