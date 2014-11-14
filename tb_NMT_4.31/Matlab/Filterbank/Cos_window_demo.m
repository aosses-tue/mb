% Cos_window_demo: script showing different windows in FFT Filterbank.

p = [];
p.audio_sample_rate	= 16000;	% Hz
p.analysis_rate		=   900;	% Hz
p.base_level		=   4/256;
p.sat_level			= 150/256;

x = [];
x = Gen_pitch_rank_tones(p, x);

window_types = { 'Hann', 'Hamming', 'Blackman', 'BlackmanHarris3', 'BlackmanHarris4'};

for n = 1:length(window_types)

	pw = p;
	pw.window = Cos_window(128, window_types{n});
	
	pw = Append_process(pw, 'FFT_filterbank_proc'); 
	pw = Append_process(pw, 'Power_sum_envelope_proc'); 
	pw = Append_process(pw, 'LGF_proc'); 

	Plot_pitch_rank_tone_response(pw, x, 1:5, window_types{n});	
end

