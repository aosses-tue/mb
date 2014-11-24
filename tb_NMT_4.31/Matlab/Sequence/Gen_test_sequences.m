% Gen_test_sequences: Script to generate example sequences for testing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODE = Implant_modes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MP1 sequence:

q_mp1.electrodes		= [  2;   4];
q_mp1.modes				= MODE.MP1;
q_mp1.current_levels	= [102; 104];
q_mp1.phase_widths		= 25.0;	% constant phase width
q_mp1.phase_gaps		=  8.0; % constant phase gap
q_mp1.periods			= 1000;	% constant period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MP2 sequence:

q_mp2.electrodes		= [ 12;  14];
q_mp2.modes				= MODE.MP2;
q_mp2.current_levels	= [202; 204];
q_mp2.phase_widths		= 25.0;	% constant phase width
q_mp2.phase_gaps		=  8.0; % constant phase gap
q_mp2.periods			= 1000;	% constant period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MP1+2 sequence:

q_mp12.electrodes		= [  3;  13; 0];
q_mp12.modes			= MODE.MP1+2;
q_mp12.current_levels	= [ 33;  43; 0];
q_mp12.phase_widths		= 25.0;	% constant phase width
q_mp12.phase_gaps		=  8.0; % constant phase gap
q_mp12.periods			= 1000;	% constant period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MP1+2 sequence with one pulse:

q1_mp12.electrodes		=  10;
q1_mp12.modes			= MODE.MP1+2;
q1_mp12.current_levels	= 100;
q1_mp12.phase_widths	= 25.0;	% constant phase width
q1_mp12.phase_gaps		=  8.0; % constant phase gap
q1_mp12.periods			= 69.4;	% constant period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MP1+2 sequence with idle pulse only:

qidle_mp12.electrodes		= 0;
qidle_mp12.modes			= MODE.MP1+2;
qidle_mp12.current_levels	= 0;
qidle_mp12.phase_widths		= 25.0;	% constant phase width
qidle_mp12.phase_gaps		=  8.0; % constant phase gap
qidle_mp12.periods			= 69.4;	% constant period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Long, fast sequence, scanning across electrodes:

e_scan	= [(1:10), (1:11), (1:12), (1:13)];	% row vector

cl_scan = [ ...
(0:255), ...
(0:2:255), (1:2:255), ...
(0:4:255), (1:4:255), (2:4:255), (3:4:255), ...
(0:8:255), (1:8:255), (2:8:255), (3:8:255), (4:8:255), (5:8:255), (6:8:255), (7:8:255)]';	% col vector

num_frames = length(cl_scan) * length(e_scan);
e	= repmat(e_scan,  length(cl_scan), 1);

qscan.electrodes		= e(:);
qscan.modes				= MODE.MP1+2;
qscan.current_levels	= repmat(cl_scan, length(e_scan), 1); 
qscan.phase_widths		= 25.0;	% constant phase width
qscan.phase_gaps		=  8.0; % constant phase gap
qscan.periods			= 69.4;	% constant period

