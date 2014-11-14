function y = IIR_filterbank_proc(p, x)

% IIR_filterbank_proc: IIR filterbank.
% y = IIR_filterbank_proc(p, x)

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
	p = Ensure_field(p, 'audio_sample_rate', 16000);
	
	p.analysis_rate = p.audio_sample_rate;
	p.sample_rate   = p.audio_sample_rate;
	
	p = Ensure_field(p, 'filter_order', 4);
	
	if ~isfield(p, 'crossover_freqs')
		p = Ensure_field(p, 'table_num', 9);
		p.crossover_freqs = ESPrit_crossover_frequencies(p.table_num, p.audio_sample_rate/2);
	end
	p.char_freqs = 0.5 * (p.crossover_freqs(1:end-1) + p.crossover_freqs(2:end));

	% Calculate derived parameters:
	
	p.num_bands = length(p.crossover_freqs) - 1;
	for n = 1:p.num_bands
		[p.numer{n}, p.denom{n}] = butter(p.filter_order ...
			,[p.crossover_freqs(n), p.crossover_freqs(n+1)]/(p.audio_sample_rate/2) ...
			,'bandpass'...
			);
	end
	y = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for n = 1:p.num_bands
		y(n,:) = filter(p.numer{n}, p.denom{n}, x)';
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
