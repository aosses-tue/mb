function result = Collate_into_sequence_test

% Collate_into_sequence_test: Test of Collate_into_sequence_proc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_chans = 6;
[m1, m2] = ndgrid(1:num_chans, 10:10:30);
env = m1 + m2;
idle = isprime(env);		% Idle for magnitudes that are prime.
skip = (rem(env, 3) == 0);	% Skip magnitudes that are divisible by 3.

% If we specify both analysis_rate and channel_stim_rate,
% then channel_stim_rate will not be quantised.
p = [];
p.num_bands				= num_chans;
p.analysis_rate			= 2400;
p.channel_stim_rate		= 10000/num_chans;

p.channel_order_type	= 'apex-to-base';
pa = Collate_into_sequence_proc(p);

p.channel_order_type	= 'base-to-apex';
pb = Collate_into_sequence_proc(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test cases:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1: No skip, no idle (e.g. CIS with continuous T stim)

qa = Collate_into_sequence_proc(pa, env);
qb = Collate_into_sequence_proc(pb, env);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. No skip, with Idle (e.g. CIS without continuous T stim)

env_idle = env;
env_idle(idle) = -1;

qai  = Collate_into_sequence_proc(pa, env_idle);
qbi  = Collate_into_sequence_proc(pb, env_idle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Skip, no idle (e.g. ACE, with base level 0)

env_skip = env;
env_skip(skip) = NaN;

qas  = Collate_into_sequence_proc(pa, env_skip);
qbs  = Collate_into_sequence_proc(pb, env_skip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Skip and Idle (e.g. ACE)

env_skip_idle = env;
env_skip_idle(idle) = -1;
env_skip_idle(skip) = NaN;

qasi  = Collate_into_sequence_proc(pa, env_skip_idle);
qbsi  = Collate_into_sequence_proc(pb, env_skip_idle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose
	env
	Disp_sequence(qa);
	Disp_sequence(qb);
	
	env_idle
	Disp_sequence(qai);
	Disp_sequence(qbi);
	
	env_skip
	Disp_sequence(qas);
	Disp_sequence(qbs);
	
	env_skip_idle
	Disp_sequence(qasi);
	Disp_sequence(qbsi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apex-to-base:

a_ = [
	1	11	-1
	2	12	12
	3	13	-1
	4	14	14
	5	15	15
	6	16	16
	1	21	21
	2	22	22
	3	23	-1
	4	24	24
	5	25	25
	6	26	26
	1	31	-1
	2	32	32
	3	33	33
	4	34	34
	5	35	35
	6	36	36
];

qa_.channels    = a_(:,1);
qa_.magnitudes  = a_(:,2);
qa_.periods		= 100;

qai_.channels   = a_(:,1);
qai_.magnitudes = a_(:,3);
qai_.periods	= 100;

% Copy a_ but comment out (skip) magnitudes that are divisible by 3
as_ = [
	1	11	-1
%	2	12	12
	3	13	-1
	4	14	14
%	5	15	15
	6	16	16
%	1	21	21
	2	22	22
	3	23	-1
%	4	24	24
	5	25	25
	6	26	26
	1	31	-1
	2	32	32
%	3	33	33
	4	34	34
	5	35	35
%	6	36	36
];

qas_.channels    = as_(:,1);
qas_.magnitudes  = as_(:,2);
qas_.periods     = 100;

qasi_.channels   = as_(:,1);
qasi_.magnitudes = as_(:,3);
qasi_.periods    = 100;

Tester(qa,		qa_);
Tester(qai,		qai_);
Tester(qas,		qas_);
Tester(qasi,	qasi_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% base-to-apex:

b_ = [
	6	16	16
	5	15	15
	4	14	14
	3	13	-1
	2	12	12
	1	11	-1
	6	26	26
	5	25	25
	4	24	24
	3	23	-1
	2	22	22
	1	21	21
	6	36	36
	5	35	35
	4	34	34
	3	33	33
	2	32	32
	1	31	-1
];
qb_.channels    = b_(:,1);
qb_.magnitudes  = b_(:,2);
qb_.periods     = 100;

qbi_.channels   = b_(:,1);
qbi_.magnitudes = b_(:,3);
qbi_.periods    = 100;

% Copy b_ but comment out (skip) magnitudes that are divisible by 3
bs_ = [
	6	16	16
%	5	15	15
	4	14	14
	3	13	-1
%	2	12	12
	1	11	-1
	6	26	26
	5	25	25
%	4	24	24
	3	23	-1
	2	22	22
%	1	21	21
%	6	36	36
	5	35	35
	4	34	34
%	3	33	33
	2	32	32
	1	31	-1
];

qbs_.channels    = bs_(:,1);
qbs_.magnitudes  = bs_(:,2);
qbs_.periods     = 100;

qbsi_.channels   = bs_(:,1);
qbsi_.magnitudes = bs_(:,3);
qbsi_.periods    = 100;

Tester(qb,		qb_);
Tester(qbi,		qbi_);
Tester(qbs,		qbs_);
Tester(qbsi,	qbsi_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result
