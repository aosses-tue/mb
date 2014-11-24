function result = Plot_sequence_test

% Plot_sequence_test: Test of Plot_sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test titles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Single seq, single title
title_str = Init_title_strings('Q', 1);
Tester(title_str, {'Q'});

% Multiple seqs, one title
title_str = Init_title_strings('A', 2);
Tester(title_str, {'A (1)','A (2)'});

% Multiple seqs, multiple titles
title_str = Init_title_strings({'A','B'}, 2);
Tester(title_str, {'A (1)','B (2)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test sequence plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All fields length N:
q1 = [];
q1.channels = (1:4)';
q1.magnitudes = 20 * q1.channels;
q1.periods = (100:50:250)';  

% Constant magnitude, period
q2 = [];
q2.channels = (1:6)';
q2.magnitudes = 50;
q2.periods = 100;

% Constant channel, period
q3 = [];
q3.channels = 8;
q3.magnitudes = (10:10:60)';
q3.periods = 100;

Plot_sequence({q1, q2, q3});

x = axis;
Tester(x(2), 700);			% maximum sequence duration (q1)
Tester(x(3:4), [0.75 9]);	% maximum channels + 1

u = Get_figure_data;
Tester(u.seqs,			{q1, q2, q3});
Tester(u.num_seqs,		3);
Tester(u.title_str,		{'Sequence (1)'  'Sequence (2)'  'Sequence (3)'});
Tester(u.min_channel,	1);
Tester(u.max_channel,	8);
Tester(u.max_mag,		80);
Tester(u.max_time,		700);

Tester(get(u.h_figure, 'Name'), 'Sequence (1)');

% Key presses:

Plot_sequence('KeyPress', ']');	% Switch to next sequence.
Tester(get(u.h_figure, 'Name'), 'Sequence (2)');

Plot_sequence('KeyPress', ']');	% Switch to next sequence.
Tester(get(u.h_figure, 'Name'), 'Sequence (3)');

Plot_sequence('KeyPress', ']');	% Wrap-around.
Tester(get(u.h_figure, 'Name'), 'Sequence (1)');

Plot_sequence('KeyPress', '[');	% Wrap-around.
Tester(get(u.h_figure, 'Name'), 'Sequence (3)');

Plot_sequence('KeyPress', '2');	% direct selection.
Tester(get(u.h_figure, 'Name'), 'Sequence (2)');

Plot_sequence('KeyPress', '0');	% Wrap-around.
Tester(get(u.h_figure, 'Name'), 'Sequence (3)');

if ~verbose
	close(gcf);
end

% Constant channel, period; with some idle pulses
q4 = q3;
q4.magnitudes(3) = -1;
Plot_sequence(q4, 'idle');
u = Get_figure_data;
Tester(get(u.h_figure, 'Name'), 'idle');

if ~verbose
	close(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result




