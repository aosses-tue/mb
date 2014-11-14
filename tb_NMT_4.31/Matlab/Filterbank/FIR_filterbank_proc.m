function y = FIR_filterbank_proc(p, x)

% FIR_filterbank_proc: FIR filterbank.
% y = FIR_filterbank_proc(p, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	y = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set default values for parameters that are absent:
	p = Ensure_field(p, 'audio_sample_rate',16000);
	
	p.analysis_rate = p.audio_sample_rate;
	p.sample_rate   = p.audio_sample_rate;
	
	p = Ensure_field(p, 'filter_order',      128);
	
	p.crossover_freqs = ESPrit_crossover_frequencies;
	p.crossover_freqs = p.crossover_freqs(3:21);
	p.char_freqs = 0.5 * (p.crossover_freqs(1:end-1) + p.crossover_freqs(2:end));

	% Calculate derived parameters:
	
	p.num_bands = length(p.crossover_freqs) - 1;
	for n = 1:p.num_bands
		p.numer{n} = fir1(p.filter_order ...
			,[p.crossover_freqs(n), p.crossover_freqs(n+1)]/(p.audio_sample_rate/2) ...
			,'bandpass'...
			);
	end
	
	y = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Convert to a row vector if necessary:
	[r,c] = size(x);
	if c == 1
		x = x.';
	end
	
	for n = 1:p.num_bands
		y(n,:) = filter(p.numer{n}, 1, x);
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
