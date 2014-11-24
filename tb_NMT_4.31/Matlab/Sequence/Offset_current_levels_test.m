function result = Offset_current_levels_test

% Offset_current_levels_test: Test of Offset_current_levels_proc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODE = Implant_modes;

q1.electrodes		= [  2;  0;  3;  1;  0;  4];
q1.modes			= MODE.MP1+2;
q1.current_levels	= [100;  0;  0;250;  0; 10];
q1.phase_widths		=  25.0;
q1.phase_gaps		=   8.0;
q1.periods			= 100.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% positive offset
% Idle pulses (defined as current_level == 0) do not change

p20.current_level_offset = 20;
p20 = Offset_current_levels_proc(p20);

r20 = Offset_current_levels_proc(p20, q1);
r20_ = q1;
r20_.current_levels	= [120;  0;  0;255;  0; 30];
Tester(r20, r20_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% negative offset

n20.current_level_offset = -20;
n20 = Offset_current_levels_proc(n20);

g20 = Offset_current_levels_proc(n20, q1);
g20_ = q1;
g20_.current_levels	= [ 80;  0;  0;230;  0;  0];
Tester(g20, g20_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result
