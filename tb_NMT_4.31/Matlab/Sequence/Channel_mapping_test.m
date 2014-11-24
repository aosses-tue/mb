function result = Channel_mapping_test;

% Channel_mapping_test: Test of Channel_mapping_proc, Channel_mapping_inv_proc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson, Herbert Mauch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODE = Implant_modes;

p1 = [];
p1.num_bands		= 22;
p1.electrodes		= (p1.num_bands:-1:1)';
p1.modes			= MODE.MP1+2;
p1.phase_width		= 25.0;		% microseconds
p1.phase_gap		=  8.0;		% microseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters chosen so that current_level = magnitude

p1.threshold_levels	= repmat(  0, p1.num_bands, 1);
p1.comfort_levels	= repmat(255, p1.num_bands, 1);
p1.full_scale		= 255;

p1 = Channel_mapping_proc(p1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep sequence has each channel in order:

N = 22;
n = (1:N)';
c0.channels			= n;
c0.magnitudes		= n * 10;
c0.periods			= repmat(100, N, 1);

q0 = Channel_mapping_proc(p1, c0);

if verbose > 1
	disp('c0');
	Disp_sequence(c0);
	disp('q0');
	Disp_sequence(q0);
end;

Tester(Get_num_pulses(c0),	Get_num_pulses(q0));
Tester(q0.electrodes,		p1.electrodes);
Tester(q0.modes,			p1.modes);
Tester(q0.current_levels,	c0.magnitudes);
Tester(q0.phase_widths,		p1.phase_width);
Tester(q0.phase_gaps,		p1.phase_gap);
Tester(q0.periods,			c0.periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse mapping:

p1_inv = Channel_mapping_inv_proc(p1);
c0_inv = Channel_mapping_inv_proc(p1_inv, q0);

Tester(c0_inv, c0);

if verbose > 1
	disp('c0_inv');
	Disp_sequence(c0_inv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence has idle pulse, then each channel in order:

c1.channels			= [ 1; c0.channels];
c1.magnitudes		= [-1; c0.magnitudes];
c1.periods			= repmat(100, N+1, 1);

q1 = Channel_mapping_proc(p1, c1);

if verbose > 1
	disp('c1');
	Disp_sequence(c1);
	disp('q1');
	Disp_sequence(q1);
end;

Tester(Get_num_pulses(c1),	Get_num_pulses(q1));
Tester(q1.electrodes,		[ 0; p1.electrodes]);
Tester(q1.modes,			p1.modes);
Tester(q1.current_levels,	[ 0; c1.magnitudes(2:end)]);
Tester(q1.phase_widths,		p1.phase_width);
Tester(q1.phase_gaps,		p1.phase_gap);
Tester(q1.periods,			c1.periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse mapping:

c1_inv = Channel_mapping_inv_proc(p1_inv, q1);

% There are several ways a (chan, mag) pair can represent an idle pulse,
% but the inverse mapping always uses (0, 0).

Tester(c1_inv.channels,		[0; c0.channels]);
Tester(c1_inv.magnitudes,	[0; c0.magnitudes]);
Tester(c1_inv.periods,		c1.periods);

if verbose > 1
	disp('c1_inv');
	Disp_sequence(c1_inv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Idle pulses:

%  Sequences with constant channel, all idle pulses:

c_idle.channels		= 3;
c_idle.magnitudes	= repmat(-1, 4, 1);
c_idle.periods		= 100;

% mix of idle and normal pulses:

c_mix.channels		= 3;
c_mix.magnitudes	= [10; 20; -1; 30; -1; -1];
c_mix.periods		= 100;

% Idle pulses without specifying channel:

c_null.channels		= [ 1;  0;  2;  0;  3;  0];
c_null.magnitudes	= [10; 20; -1; 10; 30; -1];
c_null.periods		= 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special idle:

Tester(p1.special_idle, 1);

q_special_idle = Channel_mapping_proc(p1, c_idle);
Tester(q_special_idle.electrodes,		zeros(4,1));
Tester(q_special_idle.modes,			MODE.MP1+2);
Tester(q_special_idle.current_levels,	zeros(4,1));

q_special_mix = Channel_mapping_proc(p1, c_mix);
Tester(q_special_mix.electrodes,		[20; 20;  0; 20;  0;  0]);
Tester(q_special_mix.modes,				MODE.MP1+2);
Tester(q_special_mix.current_levels,	[10; 20;  0; 30;  0;  0]);

q_special_null = Channel_mapping_proc(p1, c_null);
Tester(q_special_null.electrodes,		[22;  0;  0;  0; 20;  0]);
Tester(q_special_null.modes,			MODE.MP1+2);
Tester(q_special_null.current_levels,	[10;  0;  0;  0; 30;  0]);

if verbose > 1
	disp('q_special_idle');
	Disp_sequence(q_special_idle);
	disp('q_special_mix');
	Disp_sequence(q_special_mix);
	disp('q_special_null');
	Disp_sequence(q_special_null);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No special idles:

pn = p1;
pn.special_idle = 0;
q_normal_idle = Channel_mapping_proc(pn, c_idle);

Tester(q_normal_idle.electrodes,		20);
Tester(q_normal_idle.modes,				MODE.MP1+2);
Tester(q_normal_idle.current_levels,	zeros(4,1));

if verbose > 1
	disp('q_normal_idle');
	Disp_sequence(q_normal_idle);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters chosen so that current_level = threshold_level

p2 = p1;
p2.threshold_levels	= 5 * (1:p2.num_bands)';
p2.volume			= 0;

p2 = Channel_mapping_proc(p2);
q2 = Channel_mapping_proc(p2, c1);

Tester(Get_num_pulses(c1),	Get_num_pulses(q2));
Tester(q2.electrodes,		[ 0; p2.electrodes]);
Tester(q2.modes,			p2.modes);
Tester(q2.current_levels,	[ 0; p2.threshold_levels]);
Tester(q2.phase_widths,		p2.phase_width);
Tester(q2.phase_gaps,		p2.phase_gap);
Tester(q2.periods,			c1.periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c3.channels		= [4;3;5;1;2];
c3.magnitudes	= linspace(0,100,5)';
c3.periods		= (100:100:500)';

p3.num_bands		= 5;
p3.electrodes		= [20;16;12;7;2];
p3.modes			= MODE.MP1+2;
p3.threshold_levels	= repmat(100, p3.num_bands, 1);
p3.comfort_levels	= repmat(200, p3.num_bands, 1);
p3.full_scale		= 100;

p3.volume = 100;
p3  = Channel_mapping_proc(p3);
q3a = Channel_mapping_proc(p3, c3);

Tester(q3a.electrodes,	    [7;12;2;20;16]);
Tester(q3a.current_levels,	round(linspace(100,200,5))');
Tester(q3a.periods,			c3.periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse mapping:

p3_inv = Channel_mapping_inv_proc(p3);
c3_inv = Channel_mapping_inv_proc(p3_inv, q3a);

Tester(c3_inv, c3);

if verbose > 1
	disp('c3_inv');
	Disp_sequence(c3_inv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p3.volume = 50;
p3  = Channel_mapping_proc(p3);
q3b = Channel_mapping_proc(p3, c3);

Tester(q3b.current_levels,	round(linspace(100,150,5))');
Tester(q3b.periods,			c3.periods);

% +++ check magnitudes greater than full_scale
% +++ check volumes greater than 100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test bad volume value:

p3.volume = NaN;
try
	lasterr('');	% Clear any previous error
	q3NaN = Channel_mapping_proc(p3, c3);
	Check_sequence(p3, c3);	% would be called by Stream_sequence
	Tester(0);		% FAIL
catch
	Tester(1);
end;

if verbose > 1
	disp('Pass test if next line is an error message:');
	disp(lasterr);
	disp('.');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p4 = p3;
p4.volume_type = 'constant range';
p4.volume = 100;
p4 = Channel_mapping_proc(p4);
q4a = Channel_mapping_proc(p4, c3);

p4.volume = 50;
p4  = Channel_mapping_proc(p4);
q4b = Channel_mapping_proc(p4, c3);

p4.volume = 0;
p4  = Channel_mapping_proc(p4);
q4c = Channel_mapping_proc(p4, c3);

Tester(q4a.current_levels,	round(linspace(100,200,5))');
Tester(q4b.current_levels,	round(linspace( 50,150,5))');
Tester(q4c.current_levels,	round(linspace(  0,100,5))');

Tester(q4a.periods,			c3.periods);
Tester(q4b.periods,			c3.periods);
Tester(q4c.periods,			c3.periods);

p5.num_bands		= 5;
p5.electrodes		= (p5.num_bands:-1:1)';
p5.modes			= MODE.MP1+2;
p5.threshold_levels	= repmat(100, p5.num_bands, 1);
p5.comfort_levels	= repmat(150, p5.num_bands, 1);
p5.full_scale		= 100;
p5.volume_type		= 'constant range';
p5.volume			= 100;
p5 = Channel_mapping_proc(p5);
q5a = Channel_mapping_proc(p5, c3);

p5.volume			= 200;
p5  = Channel_mapping_proc(p5);
q5b = Channel_mapping_proc(p5, c3);

Tester(q5a.current_levels,	round(linspace(100,150,5))');
Tester(q5b.current_levels,	round(linspace(150,200,5))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence from map with disabled electrodes

p6 = [];
p6.electrodes		= [22,21,18,15,14,12,9,8,7,3]';
N           		= length(p6.electrodes);
p6.num_bands		= N;
p6.electrodes		= (2*p6.num_bands:-2:1)';
p6.modes			= MODE.MP1+2;
p6.phase_width		= 25.0;		% microseconds
p6.phase_gap		=  8.0;		% microseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters chosen so that current_level = magnitude

p6.threshold_levels	= repmat(  0, p6.num_bands, 1);
p6.comfort_levels	= repmat(255, p6.num_bands, 1);
p6.full_scale		= 255;

p6 = Channel_mapping_proc(p6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep sequence has each channel in order:

n = (1:N)';
c6.channels			= n;
c6.magnitudes		= n * 10;
c6.periods			= repmat(100, N, 1);

q6 = Channel_mapping_proc(p6, c6);

if verbose > 1
	disp('c6');
	Disp_sequence(c6);
	disp('q6');
	Disp_sequence(q6);
end;

Tester(Get_num_pulses(c6),	Get_num_pulses(q6));
Tester(q6.electrodes,		p6.electrodes);
Tester(q6.modes,			p6.modes);
Tester(q6.current_levels,	c6.magnitudes);
Tester(q6.phase_widths,		p6.phase_width);
Tester(q6.phase_gaps,		p6.phase_gap);
Tester(q6.periods,			c6.periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse mapping:

p6_inv = Channel_mapping_inv_proc(p6);
c6_inv = Channel_mapping_inv_proc(p6_inv, q6);

Tester(c6_inv, c6);

if verbose > 1
	disp('c6_inv');
	Disp_sequence(c6_inv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;
