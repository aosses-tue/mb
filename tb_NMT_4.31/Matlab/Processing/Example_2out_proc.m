function [q1, q2] = Example_2out_proc(p, x)

% Example_2out_proc: Example processing function with two outputs.
% Used in Process_chain_test to test Process_chain_nargout.
% [q1, q2] = Example_2out_proc(p, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0  % Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q1 = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1  % Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set default values for parameters that are absent:
    p = Ensure_field(p, 'h1', 3);
    p = Ensure_field(p, 'h2', 5);
    
	% Calculate derived parameters:
	p.h3 = p.h1 * p.h2;

    q1 = p;  % Return parameters.    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2  % Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q1 = x + p.h3;
    q2 = x + p.h2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
