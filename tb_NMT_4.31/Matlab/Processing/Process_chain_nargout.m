function y = Process_chain_nargout(p, x)

% Process_chain_nargout: Process a signal according to a parameter struct, return all intermediate signals.
% Includes all outputs of any functions that return multiple outputs.
%
% y = Process_chain_nargout(p, x)
%
% p:           Parameter struct.
% p.processes: Functions to call.
% x:           Input to the first process in the chain.
% y:           Cell array containing all of the outputs of each process in the chain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(p.processes)
	p = feval(p.processes{n}, p);	% Calculate parameters.
end

k = 1;
for n = 1:length(p.processes)
	f = p.processes{n};
	kend = k + nargout(f) - 1;
	[y{k:kend, 1}] = feval(f, p, x);	% Append multiple outputs to cell array.
	x = y{k};							% First output becomes next input.
	k = kend + 1;
end
