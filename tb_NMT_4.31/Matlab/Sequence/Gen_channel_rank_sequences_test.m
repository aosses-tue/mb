function result = Gen_channel_rank_sequences_test

% Gen_channel_rank_sequences_test: Test of Gen_channel_rank_sequences.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No idle pulses:

p = [];
% Choose rates so that periods are exact divisors:
p.channel_stim_rate = 2000;
p.num_selected      = 5;
p = Append_process(p, 'Channel_mapping_proc');
Tester(p.num_bands, 22);

xn = Gen_channel_rank_sequences(p, []);

Tester(length(xn.stimuli), 22);
Tester(xn.stimulus_duration, 0.50);
Tester(xn.gap_duration, 0.25);

q = xn.gap;
Tester(q.channels, 1);
Tester(Get_num_pulses(q), p.channel_stim_rate * xn.gap_duration);
Tester(q.magnitudes, repmat(-1, Get_num_pulses(q), 1));
Tester(q.periods, 1e6 / p.channel_stim_rate);

Tester(xn.stimulus_names{1}, 'Ch  1');
q = xn.stimuli{1};
Tester(q.channels, 1);
Tester(Get_num_pulses(q), p.channel_stim_rate * xn.stimulus_duration);
Tester(q.magnitudes, repmat(1, Get_num_pulses(q), 1));
Tester(q.periods, 1e6 / p.channel_stim_rate);

Tester(xn.stimulus_names{2}, 'Ch  2');
q = xn.stimuli{2};
Tester(q.channels, 2);
Tester(Get_num_pulses(q), p.channel_stim_rate * xn.stimulus_duration);
Tester(q.magnitudes, repmat(1, Get_num_pulses(q), 1));
Tester(q.periods, 1e6 / p.channel_stim_rate);

Tester(xn.stimulus_names{10}, 'Ch 10');
q = xn.stimuli{10};
Tester(q.channels, 10);
Tester(Get_num_pulses(q), p.channel_stim_rate * xn.stimulus_duration);
Tester(q.magnitudes, repmat(1, Get_num_pulses(q), 1));
Tester(q.periods, 1e6 / p.channel_stim_rate);

Tester(xn.stimulus_names{22}, 'Ch 22');
q = xn.stimuli{22};
Tester(q.channels, 22);
Tester(Get_num_pulses(q), p.channel_stim_rate * xn.stimulus_duration);
Tester(q.magnitudes, repmat(1, Get_num_pulses(q), 1));
Tester(q.periods, 1e6 / p.channel_stim_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insert idle pulses explicitly:

xi.insert_idles = 1;
xi = Gen_channel_rank_sequences(p, xi);

Tester(p.implant_stim_rate, 10000);
Tester(length(xi.stimuli), 22);
Tester(xi.stimulus_duration, 0.50);
Tester(xi.gap_duration, 0.25);

q = xi.gap;
Tester(q.channels, 1);
Tester(Get_num_pulses(q), p.implant_stim_rate * xi.gap_duration);
Tester(q.magnitudes, repmat(-1, Get_num_pulses(q), 1));
Tester(q.periods, 100);

Tester(xi.stimulus_names{7}, 'Ch  7');
q = xi.stimuli{7};
Tester(q.channels, 7);
Tester(Get_num_pulses(q), p.implant_stim_rate * xi.stimulus_duration);
% 2000 Hz for 0.5 seconds = 1000 stim pulses.
% interleaved with (p.num_selected - 1) = 4 idle pulses
Tester(q.magnitudes, repmat([1;-1;-1;-1;-1], 1000, 1));
Tester(q.periods, 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insert idle pulses by default:

p250 = [];
% Choose rates so that periods are exact divisors:
p250.channel_stim_rate = 250;
p250.num_selected      = 8;
p250 = Append_process(p250, 'Channel_mapping_proc');

x250 = Gen_channel_rank_sequences(p250, []);

Tester(p250.implant_stim_rate, 2000);

q = x250.gap;
Tester(q.channels, 1);
Tester(Get_num_pulses(q), 500);
Tester(q.magnitudes, repmat(-1, 500, 1));
Tester(q.periods, 500);

q = x250.stimuli{7};
Tester(q.channels, 7);
% 2000 pps for 0.5 seconds = 1000 total pulses.
Tester(Get_num_pulses(q), 1000);
% 250 pps for 0.5 seconds = 125 stim pulses.
% interleaved with (p.num_selected - 1) = 7 idle pulses
Tester(q.magnitudes, repmat([1;-1;-1;-1;-1;-1;-1;-1], 125, 1));
Tester(q.periods, 500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Omit idle pulses explicitly:

xn250.insert_idles = 0;
xn250 = Gen_channel_rank_sequences(p250, xn250);

q = xn250.gap;
Tester(q.channels, 1);
Tester(Get_num_pulses(q), 63);
Tester(q.magnitudes, repmat(-1, 63, 1));
Tester(q.periods, 4000);

q = xn250.stimuli{7};
Tester(q.channels, 7);
% 250 pps for 0.5 seconds = 125 stim pulses.
Tester(Get_num_pulses(q), 125);
Tester(q.magnitudes, repmat(1, 125, 1));
Tester(q.periods, 4000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result
