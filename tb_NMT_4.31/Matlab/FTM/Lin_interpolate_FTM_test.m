function result = Lin_interpolate_FTM_test

% Lin_interpolate_FTM_test: Test of Lin_interpolate_FTM_proc

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

m1 = Lin_interpolate_FTM_proc(1, m);

% With parameter struct:
p1.analysis_rate     = 1000;
p1.channel_stim_rate = p1.analysis_rate;
v1 = Lin_interpolate_FTM_proc(p1, m);

if verbose
	m1, v1
end

Tester(m1, m);
Tester(v1, m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double the rate:

m2 = Lin_interpolate_FTM_proc(2, m);

% With parameter struct:
p2.analysis_rate     = 1000;
p2.channel_stim_rate = 2 * p2.analysis_rate;
v2 = Lin_interpolate_FTM_proc(p2, m);

x = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80];
m2_ = [x; x+1];

if verbose
	m2, v2
end

Tester(m2, m2_);
Tester(v2, m2_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-integral multiple:

m3 = Lin_interpolate_FTM_proc(3/2, m);

% With parameter struct:
p3.analysis_rate     =  1000;
p3.channel_stim_rate = (3/2) * p3.analysis_rate;
v3 = Lin_interpolate_FTM_proc(p3, m);

% leave off 1 input:
ma  = m(:,1:(end-1));
m3a = Lin_interpolate_FTM_proc(3/2, ma);
v3a = Lin_interpolate_FTM_proc(p3,  ma);

if verbose
	m3, v3, m3a, v3a
end

x = [10, 10+20/3, 20+10/3, 30, 30+20/3, 40+10/3, 50, 50+20/3, 60+10/3, 70, 70+20/3];
m3_ = [x; x+1];

x(end) = [];
m3a_ = [x; x+1];

Tester(m3,  m3_,  1e-10);
Tester(v3,  m3_,  1e-10);
Tester(m3a, m3a_, 1e-10);
Tester(v3a, m3a_, 1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result

