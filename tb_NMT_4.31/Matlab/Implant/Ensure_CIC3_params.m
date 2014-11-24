function implant = Ensure_CIC3_params(implant)

% Ensure_CIC3_params: Ensure implant parameters are valid for CIC3.
%
% p_out = Ensure_CIC3_params(p_in)
%
% p_in:  A struct containing the implant parameters.
%          Any fields omitted will be set to default values.
% p_out: A struct containing the clinical and derived parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    implant = [];
end

implant.IC = 'CIC3';

implant = Ensure_field(implant, 'protocol', 'embedded');
implant = Ensure_field(implant, 'rf_freq',  5e6);			% Hz

implant = Embedded_protocol_params(implant);

implant.minimum_phase_width		= 24.6;
implant.default_phase_width   	= 25.0;
implant.maximum_phase_width		= 75.6;	% Maximum supported by SPrint streaming

implant.default_phase_gap     	=  8.0;

implant.CURRENT_LEVEL_MAX		= 255;

implant.CURRENT_uA_MIN			=   10.00;
implant.CURRENT_uA_MAX			= 1750.00;
