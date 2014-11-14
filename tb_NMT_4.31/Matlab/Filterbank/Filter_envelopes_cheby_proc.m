function y = Filter_envelopes_cheby_proc(p, x)

% Filter_envelopes_cheby_proc: Apply a Chebychev IIR filter to each envelope.
% y = Filter_envelopes_cheby_proc(p, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Cochlear Ltd
%   $Header: //dsp/Nucleus/_Releases/NMT_4.30/Matlab/Filterbank/Filter_envelopes_cheby_proc.m#1 $
% $DateTime: 2008/03/04 14:27:13 $
%   $Change: 86418 $
%   Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	y = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set default values for parameters that are absent:
	p = Ensure_field(p, 'analysis_rate',	16000);
	p = Ensure_field(p, 'env_order',            4);
	p = Ensure_field(p, 'env_ripple',      		0.5);
	p = Ensure_field(p, 'env_pass_band_Hz',   400);
	p = Ensure_field(p, 'env_type',    	    'low');

	% Calculate filter coefficients:
	[p.env_numer, p.env_denom] = cheby1(p.env_order, p.env_ripple,...
			p.env_pass_band_Hz/(p.analysis_rate/2), p.env_type);
	
	y = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Apply filter to each envelope:
	y = filter(p.env_numer, p.env_denom, x, [], 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
