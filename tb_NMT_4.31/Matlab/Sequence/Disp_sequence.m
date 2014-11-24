function Disp_sequence(seq)
% Disp_sequence: Displays a sequence struct as text in the command window.
% function Disp_sequence(seq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f, old_echo] = Open_log;
Open_log('echo');	% Turn echo on.
Log_sequence(seq);
Open_log(old_echo);	% Restore old echo state.
