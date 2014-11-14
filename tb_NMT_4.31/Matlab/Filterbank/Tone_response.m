function [output, response, max_gain] = Tone_response(p, freq_vec, num_in_samples)

% Tone_response: Calculate frequency response of filterbank using tones (& plot if verbose)
%
% [output, response, max_gain] = Tone_response(p, freq_vec, num_in_samples)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (0:(num_in_samples - 1))';

if length(freq_vec) == 1
	freq_vec = 0:freq_vec:(p.audio_sample_rate/2);
end

num_tones = length(freq_vec);
for k = 1:num_tones
	tone = sin(2*pi*n* freq_vec(k)/p.audio_sample_rate);
	output{k} = Process(p, tone);
	steady_output = output{k};
	response(:,k) = max(abs(steady_output), [], 2);
end

max_gain = max(response, [], 2);

if Verbose >= 1	
	disp('Max gain of each band =');
	disp(max_gain);
end
if Verbose >= 2
	Plot_freq_response(freq_vec, response, 'log');
end
if Verbose >= 3
	Plot_freq_response(freq_vec, response, 'linear');
	GUI_FTM(p, response, 'Response (FTM)');	
	GUI_FTM(p, output,   'Output');
end
