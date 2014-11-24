% Gain_demo: Demonstrates effect of varying gain on compressed filter envelopes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

speech = wavread('choice');
speech_rms = sqrt(mean(speech .* speech));

randn('state',0);
noise = 0.01 * randn(size(speech));
noise_rms = sqrt(mean(noise .* noise));

snr = 20*log10(speech_rms/noise_rms);

audio = speech + noise;

if Verbose
	speech_rms
	noise_rms
	snr
	figure;
	plot([speech;noise;audio]);
end

p = [];
p = Append_process(p, 'FFT_filterbank_proc'); 
p = Append_process(p, 'Power_sum_envelope_proc'); 
p = Append_process(p, 'LGF_proc'); 

gaudio = [];
titles = [];
envelope = [];
gain_dB = -20:4:40;
for n = 1:length(gain_dB)
	gain        = 10 ^ (gain_dB(n) / 20);
	titles{n}   = sprintf('%d dB', gain_dB(n));
	gaudio{n}   = audio * gain;
	envelope{n} = Process(p, gaudio{n});
end

GUI_FTM(p, envelope, titles);


	
	