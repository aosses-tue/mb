function result = Ensure_rate_params_test

% Ensure_rate_params_test: Test Ensure_rate_params. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantisation of analysis rate & channel stim rate:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There are four cases:

% Specify neither rate:
%	both are set equal to a default (quantised) value.

p0 = Ensure_rate_params;
Tester(p0.analysis_rate,		500);
Tester(p0.channel_stim_rate,	500);

% Specify channel_stim_rate only (most common):
%	analysis_rate is set equal to channel_stim_rate,
%	then both are quantised.

p1.channel_stim_rate	=		1005;
p1.num_selected			=		  12;
p1 = Ensure_rate_params(p1);
Tester(p1.analysis_rate,		1000);
Tester(p1.channel_stim_rate,	1000);

% Specify analysis_rate only:
%	analysis_rate is quantised.
%	channel_stim_rate is set equal to quantised analysis_rate.

p2.analysis_rate		=		1005;
p2.num_selected			=		  12;
p2 = Ensure_rate_params(p2);
Tester(p2.analysis_rate,		1000);
Tester(p2.channel_stim_rate,	1000);

% Specify both analysis_rate & channel_stim_rate:
%	analysis_rate is quantised.
%	channel_stim_rate is not changed.
%	This will (usually) require Resample_FTM_proc.

p3.analysis_rate		=		1005;
p3.channel_stim_rate	=		2400;
p3.num_selected			=		   6;
p3 = Ensure_rate_params(p3);
Tester(p3.analysis_rate,		1000);
Tester(p3.channel_stim_rate,	2400);

p4.audio_sample_rate	=		16000;
p4.analysis_rate		=		p4.audio_sample_rate;
p4.channel_stim_rate	=		1800;
p4.num_selected			=		   8;
p4 = Ensure_rate_params(p4);
Tester(p4.analysis_rate,		16000);
Tester(p4.channel_stim_rate,	1800);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check rate against implant capability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default to CIC3:
pr0.audio_sample_rate	= 	  1e6/64;
pr0.channel_stim_rate	=		2000;
pr0.num_selected		=		  12;
try
	lasterr('');	% Clear any previous error
	pr0 = Ensure_rate_params(pr0);
	Tester(0);		% FAIL if we reach here.
catch
	Tester(String_contains(lasterr, 'implant stimulation rate exceeds maximum'));
end

% Explicitly select CIC3:
pr3.audio_sample_rate	= 	  1e6/64;
pr3.implant.IC			=     'CIC3';
pr3.channel_stim_rate	=		2000;
pr3.num_selected		=		  12;
try
	lasterr('');	% Clear any previous error
	pr3 = Ensure_rate_params(pr3);
	Tester(0);		% FAIL if we reach here.
catch
	Tester(String_contains(lasterr, 'implant stimulation rate exceeds maximum'));
end

% Explicitly select CIC4:
pr4.audio_sample_rate	= 	  1e6/64;
pr4.implant.IC			=     'CIC4';
pr4.channel_stim_rate	=		pr4.audio_sample_rate/6;
pr4.num_selected		=		  12;
pr4 = Ensure_rate_params(pr4);
Tester(pr4.analysis_rate,		1e6/(6*64));
Tester(pr4.channel_stim_rate,	1e6/(6*64));
Tester(pr4.implant_stim_rate,	2 * pr4.audio_sample_rate);
Tester(pr4.period,				32);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;
