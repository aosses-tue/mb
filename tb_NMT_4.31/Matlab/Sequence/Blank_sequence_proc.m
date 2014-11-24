function q = Blank_sequence_proc(p, q)

% Blank_sequence_proc: Periodically "blank out" sections of a sequence.
% The pulses in the "blanked" section are replaced by idle pulses.
%
% q_out = Blank_sequence_proc(p, q_in)
%
% p:                Parameter struct:
% p.blank_duration:   Duration of each blank interval (seconds).
% p.blank_freq:       Number of blank intervals per second.
% q_in:             Input pulse sequence.
% q_out:            Output pulse sequence.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Tim Neal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	q = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set default values for parameters that are absent:
	p = Ensure_field(p, 'blank_duration', 0.05);
	p = Ensure_field(p, 'blank_freq', 5);

    % Calculate derived parameters:
	p.num_frames_blanked = ceil(p.blank_duration*p.implant_stim_rate);
    p.blanking_frame_period = ceil((1/p.blank_freq)*p.implant_stim_rate);
    
	q = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:p.blanking_frame_period:length(q.electrodes)
        q.electrodes(i:(i+p.num_frames_blanked-1)) = 0;
        q.current_levels(i:(i+p.num_frames_blanked-1)) = 0;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
