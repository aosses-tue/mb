function v = Reject_smallest_proc(p, v)

% Reject_smallest_proc: rejects the smallest input magnitudes.
% function v = Reject_smallest_proc(p, v)
%
% Inputs:
% p:            The number of magnitudes to reject in each column,
%               or a parameter struct with fields:
% p.num_bands:     The number of filter bands.
% p.num_selected:  The number of magnitudes to keep in each column.
% v:            Envelope FTM (magnitudes).
%
% Outputs:
% v:            An FTM (same size as input FTM),
%               with the smallest values in each column rejected,
%               by replacing them with NaN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	v = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set default values for parameters that are absent:
	p = Ensure_field(p, 'num_bands',  22);
	p = Ensure_field(p, 'num_selected', min(p.num_bands, 12));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	v = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[num_bands num_time_slots] = size(v);
	
	if isstruct(p)
		num_rejected = num_bands - p.num_selected;
	else
		num_rejected = p;
	end

	if (num_rejected > 0)

		if (num_rejected > num_bands)
			error('Number to be rejected is greater than number of bands');
		end

		% If we treat the input matrix v as one long column vector,
		% then the indexes of the start of each column are:
		x0 = num_bands * [0:(num_time_slots-1)];

		for n = 1:num_rejected

			[m k] = min(v);

			% m is x row vector containing the minimum of each column (time-slot).
			% k is x row vector containing the row number (channel number) of each minimum.
			% If we treat the input matrix v as one long column vector,
			% then the indexes of the minima are:

			xk = x0 + k;

			v(xk) = NaN;	% Reject the smallest values.

		end

	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
