function cseq = Channel_mapping_inv_proc(p, pseq)

% Channel_mapping_inv_proc: Convert a pulse sequence back to a channel-magnitude sequence. 
%
% cseq = Channel_mapping_inv_proc(p, pseq)
%
% p:        client map (parameter struct).
% pseq:		pulse sequence.
% cseq:     channel-magnitude sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson, Herbert Mauch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	cseq = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = Channel_mapping_proc(p);
	    
	cseq = p;	% return parameters
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Create the fields in the order we like to show them in Disp_sequence.
	
    % Append idle info:
	electrodes			= [p.electrodes;		0];
	threshold_levels	= [p.threshold_levels;	0];
	comfort_levels		= [p.comfort_levels;	1];
    
	% Channels:
	
	% Electrode == 0 indicates idle pulses:
 	idle_pulses = (pseq.electrodes == 0);
	% Temporarily map idles to idle channel,
	% so that we can look up a temporary T & C:
	pseq.electrodes(idle_pulses) = electrodes(end);
	
    % Electrode to channel allocation
    
    for n=1:length(pseq.electrodes)
       cseq.channels(n,1) = find(electrodes == pseq.electrodes(n,1));
    end
    
	% Magnitude:

	volume			= p.volume/100;
	p_ranges		= comfort_levels - threshold_levels;

	q_t = threshold_levels(cseq.channels);
	q_r = p_ranges        (cseq.channels);
	
	switch (p.volume_type)

		case 'standard'	
			cseq.magnitudes = (pseq.current_levels - q_t) ./ (q_r .* volume);

		case 'constant range'
			cseq.magnitudes = ((pseq.current_levels - q_t) ./ q_r) + (1 - volume);

		otherwise
			error('Bad volume_type');
	end

	cseq.magnitudes = p.full_scale * cseq.magnitudes;

	cseq.channels  (idle_pulses) = 0;
	cseq.magnitudes(idle_pulses) = 0;

	cseq.periods = pseq.periods;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
