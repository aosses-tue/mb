function Window_title(title_str)
% Window_title: Sets the title of the window of the current figure.
% function Window_title(title_str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 176482 $
%    $Revision: #1 $
%    $DateTime: 2011/09/08 12:56:31 $
%      Authors: Brett Swanson
%  $Nokeywords: $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'numbertitle', 'off');
set(gcf, 'name', title_str);

