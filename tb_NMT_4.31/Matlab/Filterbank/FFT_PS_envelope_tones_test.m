function result = FFT_PS_envelope_tones_test

% FFT_PS_envelope_tones_test: Test of FFT_filterbank_proc & Power_sum_envelope_proc.
% To measure steady-state gains to pure tones,
% the length of the signal can be the same as the FFT length,
% and the analysis rate can be chosen to have no overlap between blocks,
% so the output will be a single column.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 128;

p.audio_sample_rate	= 16000;	% Hz
p.analysis_rate		= p.audio_sample_rate/N;
p.num_bands			= 22;
p = Append_process(p, 'FFT_filterbank_proc'); 
p = Append_process(p, 'Power_sum_envelope_proc'); 

n = (0:N-1)';

Fs = p.bin_freq/8;
freqs = 0: Fs: 7500;
for k = 1:length(freqs)
	tone = sin(2*pi*n* freqs(k)/p.audio_sample_rate);
	response(:,k) = Process(p, tone);
end

max_gain = max(response, [], 2);

if verbose
	p.sample_rate = 1000/Fs;
	GUI_FTM(p, response, 'Freq response');
	disp('Max gain of each band =');
	disp(max_gain);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tester(size(response), [p.num_bands, length(freqs)]);

tol = 1/65535;
Tester(max_gain, ones(p.num_bands, 1), tol);

% Compare to "master" response:
file_name = 'tone_response';
master_response = Read_FTM(file_name);
Tester(response, master_response, tol);

result = Tester;	% Report result

