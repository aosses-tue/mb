% FFT_VS_ESPrit_emulation_demo: Emulate ESPrit filterbank with FFT vector sum filterbank.
%
% Plots the calculated frequency response for each table.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_nums = 1:15;
for n = 1:length(table_nums)
	p = [];
	p.audio_sample_rate = 16000;
	p.crossover_freqs = ESPrit_crossover_frequencies(table_nums(n), p.audio_sample_rate/2);
	p = Append_process(p, 'FFT_filterbank_proc'); 
	p = Append_process(p, 'Vector_sum_proc');
	Plot_freq_response(p.response_freqs, p.freq_response, 'log');
	Window_title(sprintf('Table %d', table_nums(n)));
end
