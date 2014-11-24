function q = Channel_mapping_proc(p, cseq)

% Channel_mapping_proc: Map a channel-magnitude sequence to a pulse sequence.
%
% q = Channel_mapping_proc(p, cseq)
%
% p:        client map (parameter struct).
% cseq:     channel-magnitude sequence.
%
% q:		pulse sequence.

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

	q = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	MODE = Implant_modes;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Defaults:
			
	p = Ensure_implant_params(p);
	p = Ensure_field(p,'electrodes',		(22:-1:1)');	
	p = Ensure_field(p,'modes',				MODE.MP1+2);	

	if length(p.modes) > 1
		error('Only constant mode is supported');
	end	
	switch (p.modes)

		case {MODE.MP1, MODE.MP2}
			p = Ensure_field(p,'special_idle',	0);
			
		case MODE.MP1+2
			p = Ensure_field(p,'special_idle',	1);
		
		otherwise
			error('Unsupported mode');
	end

	p = Ensure_field(p,'num_bands',			length(p.electrodes));	

	p = Ensure_field(p,'threshold_levels',	repmat(  1, p.num_bands, 1));
	p = Ensure_field(p,'comfort_levels',	repmat(100, p.num_bands, 1));
	
	if ~all(p.comfort_levels >= p.threshold_levels)
		error('Threshold level exceeds comfort level');
	end
		
	p = Ensure_field(p,'phase_width',		p.implant.default_phase_width);	% microseconds
	p = Ensure_field(p,'phase_gap',			p.implant.default_phase_gap);	% microseconds
	
	p = Ensure_field(p,'full_scale',  		   1.0);
	p = Ensure_field(p,'volume',			   100);
	p = Ensure_field(p,'volume_type',		'standard');
	
	p = Ensure_rate_params(p);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	q = p;	% Return parameters.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~all(cseq.channels >= 0) | ~all(cseq.channels <= p.num_bands)
		error('Channel number out of range');
	end

	% Append idle info:
	
	electrodes			= [p.electrodes;		0];
	threshold_levels	= [p.threshold_levels;	0];
	comfort_levels		= [p.comfort_levels;	0];

	idle_pulses = (cseq.channels == 0);
	cseq.channels(idle_pulses) = length(electrodes);

	% Create the fields in the order we like to show them in Disp_sequence:

	q.electrodes	= electrodes(cseq.channels);
	q.modes			= p.modes;					% Constant mode.
	
	% Current level:
	
	volume			= p.volume/100;
	ranges			= comfort_levels - threshold_levels;
	
	q_magnitudes	= cseq.magnitudes / p.full_scale;
	q_magnitudes	= min(q_magnitudes, 1.0);

	q_t = threshold_levels(cseq.channels);
	q_r = ranges(cseq.channels);
		
	switch (p.volume_type)

		case 'standard'	
			q.current_levels = round(q_t + q_r .* volume .* q_magnitudes);

		case 'constant range'
			q.current_levels = round(q_t + q_r .* (q_magnitudes + volume - 1));

		otherwise
			error('Bad volume_type');
	end
	
	% Idle pulses are marked by magnitudes less than zero.
	q_is_idle = (q_magnitudes < 0);		% logical vector
	% The current levels calculated above do not apply for idle pulses,
	% set idle pulses current level to zero:
	q.current_levels(q_is_idle)	= 0;
	
	% If special_idle, then also set the idle pulses electrode to zero:
	if (p.special_idle) & any(q_is_idle)
		if length(q.electrodes) == 1	% replicate constant electrode
			q.electrodes = repmat(q.electrodes, length(q_is_idle), 1);
		end
		q.electrodes(q_is_idle) = 0;
	end
	
	q.phase_widths		= p.phase_width;		% Constant phase width.
	q.phase_gaps		= p.phase_gap;			% Constant phase gap.
	
	q.periods			= cseq.periods;			% Copy periods.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
