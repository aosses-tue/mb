function result = Overlap_sequences_test

% Overlap_sequences_test: Test of Overlap_sequences.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test sequences:

u.channels		= (1:4)';
u.magnitudes	= 100 + 10 * u.channels;
u.periods		= 400;

v.channels		= (11:14)';
v.magnitudes	= 100 + 10 * v.channels;
v.periods		= 400;

w.channels		= (15:18)';
w.magnitudes	= 100 + 10 * w.channels;
w.periods		= [400; 400; 400; 200];

x.channels		= [5; 6];
x.magnitudes	= 100 + 10 * x.channels;
x.periods		= [100; 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay 0 gives collisions (period 0):
% (can use a separate function to resolve collisions)

uv0 = Overlap_sequences(u,v,0);

uv0_.channels	= [ 1; 11; 2; 12; 3; 13; 4; 14];
uv0_.magnitudes	= 100 + 10 * uv0_.channels;
uv0_.periods	= [0;400;0;400;0;400;0;400];

Tester(uv0,uv0_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay shorter than first period:

uva = Overlap_sequences(u,v,100);

uva_.channels	= [ 1; 11; 2; 12; 3; 13; 4; 14];
uva_.magnitudes = 100 + 10 * uva_.channels;
uva_.periods	= [100;300;100;300;100;300;100;400];

Tester(uva,uva_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay shorter than first period:

uvb = Overlap_sequences(u,v,200);

uvb_ = uva_;
uvb_.periods	= [200;200;200;200;200;200;200;400];

Tester(uvb,uvb_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delay longer than first period:

uvc = Overlap_sequences(u,v,500);

uvc_.channels	= [ 1; 2; 11; 3; 12; 4; 13; 14];
uvc_.magnitudes = 100 + 10 * uvc_.channels;
uvc_.periods	= [400;100;300;100;300;100;400;400];

Tester(uvc,uvc_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative delay:

uvd = Overlap_sequences(u,v,-200);
vud = Overlap_sequences(v,u, 200);

uvd_.channels	= [11; 1; 12; 2; 13; 3; 14; 4];
uvd_.magnitudes = 100 + 10 * uvd_.channels;
uvd_.periods	= [200;200;200;200;200;200;200;400];

Tester(uvd,uvd_);
Tester(vud,uvd_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative delay so that 2 q2 pulses precede first q1 pulse:

uve = Overlap_sequences(u,v,-600);
vue = Overlap_sequences(v,u, 600);

uve_.channels	= [11; 12; 1; 13; 2; 14; 3; 4];
uve_.magnitudes = 100 + 10 * uve_.channels;
uve_.periods	= [400;200;200;200;200;200;400;400];

Tester(uve,uve_);
Tester(vue,uve_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w has short last period, so overlap has uniform period:

uwa = Overlap_sequences(u,w,200);

uwa_.channels	= [1; 15; 2; 16; 3; 17; 4; 18];
uwa_.magnitudes = 100 + 10 * uwa_.channels;
uwa_.periods	= [200;200;200;200;200;200;200;200];

Tester(uwa,uwa_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different length sequences:

wxa = Overlap_sequences(w,x,100);

wxa_.channels	= [15; 5; 6; 16; 17; 18];
wxa_.magnitudes = 100 + 10 * wxa_.channels;
wxa_.periods	= [100;100;200;400;400;200];

Tester(wxa,wxa_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose
	Plot_sequence({u,v,uv0,uva,uvb,uvc,uvd,uve},{'u','v','uv0','uva','uvb','uvc','uvd','uve'});
	Plot_sequence({u,w,uwa},{'u','w','uwa'});
	Plot_sequence({w,x,wxa},{'w','x','wxa'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result
