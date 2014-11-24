% ACE_chirp_demo: Electrodogram of ACE with chirp input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [];
p.audio_sample_rate	= 16000;
p.analysis_rate		=  1800;
p.num_selected		= 8;
p.num_bands			= 22;
p.Q					= 30;
p.base_level		= 8/256;
p = ACE_map(p);

% Chirp
duration = 2;
freq1 =    0;
freq2 = 8000;
rise_time = 0.1;
t  = 0:(1/p.audio_sample_rate):duration;
ch = chirp(t, freq1, duration, freq2)';

n_rise = round(rise_time * p.audio_sample_rate);	% Number of time samples.
n_flat = length(t) - 2 * n_rise;
env = Gen_smooth_rise_fall(n_rise, n_flat);
audio = ch .* env;

q = Process(p, audio);

Plot_sequence(q, 'Chirp', 1:22);

