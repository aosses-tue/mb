% GUI_FTM_demo: Demonstrate GUI to display a Frequency-Time-indexed Matrix (FTM).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F0 = 125;
p.sample_rate = 16000;
N  = 1024;							% Number of time samples.
num_tones = 6;
p.char_freqs = (1:num_tones) * F0;

n = -(N/2):(N/2)-1;					% Time index, row vector.
g = exp(-n .* n * 8 / (N*N));		% Gaussian envelope
tones = zeros(num_tones, N);
for k = 1:num_tones,
	tones(k,:) = g .* exp(j*(2*pi*n * p.char_freqs(k) / p.sample_rate));
end

GUI_FTM(p, tones, 'tones');
GUI_FTM(p, tones(:, 32:65), 'tones');