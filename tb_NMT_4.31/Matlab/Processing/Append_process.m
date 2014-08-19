function p = Append_process(p, f, branch)

% Append_process: Append a processing function to a chain.
% This function sets up the parameter struct for Process.
%
% p_out = Append_process(p_in, f, branch)
%
% p_in:     Input process parameter struct.
% f:        Processing function.
% branch:   Optional name of process branch (else uses main processes chain).
% p_out:    Output process parameter struct.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
	branch = 'processes';
end
p = Ensure_field(p, branch, {});
p.(branch){end + 1, 1} = f;				% Append function to process list.	
p = feval(f, p);						% Set up parameters.	

