function y = Process_chain(p, x)
% function y = Process_chain(p, x)
%
% Process a signal according to a parameter struct, return intermediate signals.
%
% y = Process_chain(p, x)
%
% p:           Parameter struct.
% p.processes: Functions to call.
% x:           Input to the first process in the chain.
% y:           Cell array containing the outputs of each process in the chain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(p.processes)
	p = feval(p.processes{n}, p);		% Calculate parameters.
end

in = x;
for n = 1:length(p.processes)
	y{n, 1} = feval(p.processes{n}, p, in);
	in = y{n};
end
