function u = Get_figure_data

% Get_figure_data: Get data of current callback figure or else current figure.
% An error occurs if there is no figure.
% A GUI that uses this function can have its callbacks
% called from the command line or another function, allowing 
% automated testing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_handle = Get_figure;
u = get(fig_handle, 'UserData');

