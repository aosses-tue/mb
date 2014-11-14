function result = Resample_FTM_test

% Resample_FTM_test: Test of Resample_FTM_proc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 10:10:80;
m = [x; x+1];
if verbose
	m
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equal numbers of input & output samples:

m1 = Resample_FTM_proc(1, m);

% With parameter struct:
p1.analysis_rate     = 1000;
p1.channel_stim_rate = p1.analysis_rate;
v1 = Resample_FTM_proc(p1, m);

if verbose
	m1, v1
end

Tester(m1, m);
Tester(v1, m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Halve the number of samples:

mh = Resample_FTM_proc(0.5, m);

% With parameter struct:
ph.analysis_rate     = 1000;
ph.channel_stim_rate = 0.5 * ph.analysis_rate;
vh = Resample_FTM_proc(ph, m);

x = 10:20:80;
mh_ = [x; x+1];

if verbose
	mh, vh
end

Tester(mh, mh_);
Tester(vh, mh_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double the number of samples:

m2 = Resample_FTM_proc(2, m);

% With parameter struct:
p2.analysis_rate     = 1000;
p2.channel_stim_rate = 2 * p2.analysis_rate;
v2 = Resample_FTM_proc(p2, m);

x = [10 10 20 20 30 30 40 40 50 50 60 60 70 70 80 80];
m2_ = [x; x+1];

if verbose
	m2, v2
end

Tester(m2, m2_);
Tester(v2, m2_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-integral multiple:

m3 = Resample_FTM_proc(3/2, m);

% With parameter struct:
p3.analysis_rate     =  1000;
p3.channel_stim_rate = (3/2) * p3.analysis_rate;
v3 = Resample_FTM_proc(p3, m);

% leave off 1 input:
ma  = m(:,1:(end-1));
m3a = Resample_FTM_proc(3/2, ma);
v3a = Resample_FTM_proc(p3,  ma);

if verbose
	m3, v3, m3a, v3a
end

x = [10 10 20 30 30 40 50 50 60 70 70 80];
m3_ = [x; x+1];

x(end) = [];
m3a_ = [x; x+1];

Tester(m3,  m3_);
Tester(v3,  m3_);
Tester(m3a, m3a_);
Tester(v3a, m3a_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

