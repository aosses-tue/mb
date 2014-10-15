function varargout = Process(p, x, recalc)

% Process: Process a signal according to a parameter struct. 
% function [x, p] = Process(p, x)
%
% Inputs:
% p:           Parameter struct.
% p.processes: Cell array containing a list of functions to call.
% x:           Input to the first function in the process chain.
% recalc:      If true, then parameters are recalculated before
%                processing the input signal.
%                Defaults to true if the argument is omitted.
%
% Outputs:
% x:           Output from the last function in the process chain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

case 0
		error('First arg must be a parameter struct.');	
case 1
		p = Process_parameters(p);
		varargout = {p};	
case 2
		p = Process_parameters(p);
		y = Process_signal(p, x);
		varargout = {y, p};	
case 3
		if (recalc)
			p = Process_parameters(p);
		end
		y = Process_signal(p, x);
		varargout = {y, p};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = Process_parameters(p)

	for n = 1:length(p.processes)
		p = feval(p.processes{n}, p);		% Calculate parameters.
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = Process_signal(p, x)

	if ~isfield(p, 'input')
	
		for n = 1:length(p.processes)
			x = feval(p.processes{n}, p, x);	% Perform processing.
		end
		y = x;
		
	else
	
		y = [];
		while 1
			u = feval(p.input, p, x);
			if isempty(u)
				break;
			end

			for n = 1:length(p.processes)
				u = feval(p.processes{n}, p, u);    % Perform processing.
			end

			y = feval(p.output, y, u);
		end
		
	end