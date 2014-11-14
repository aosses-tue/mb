function v = Resample_FTM_proc(p, u)

% Resample_FTM_proc: resamples a Frequency-Time Matrix to the channel stimulation rate.
% The resampling is done without interpolation, so each stimulation frame is based on
% the results of the last performed FFT
%
% function v = Resample_FTM_proc(p, u)
%
% Inputs:
% p:                  Parameter struct with the following fields:
% p.analysis_rate:      Analysis_rate of the filterbank.
% p.channel_stim_rate:  Channel stimulation rate.
% u:                  FTM sampled at the analysis rate.
%
% Outputs:
% v:                  FTM sampled at the channel stimulation rate.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Bas van Dijk, Brett Swanson
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
	p = Ensure_field(p, 'analysis_rate',     500);
	p = Ensure_field(p, 'channel_stim_rate', p.analysis_rate);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	p.sample_rate = p.channel_stim_rate;	% for GUI_FTM

	v = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if isstruct(p)
		output_factor = p.channel_stim_rate / p.analysis_rate;
	else
		output_factor = p;
	end
	
	if (output_factor == 1)
		v = u;
		return;
	end
	
	num_samples_in  = size(u, 2);
	num_samples_out = round(num_samples_in * output_factor);
	sample_points   = (0:num_samples_out-1) / output_factor;
	sample_indices  = round(sample_points + 0.5);

	v = u(:, sample_indices);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
