function implant = Ensure_CIC4_params(implant)

% Ensure_CIC4_params: Ensure implant parameters are valid for CIC4.
%
% p_out = Ensure_CIC4_params(p_in)
%
% p_in:  A struct containing the implant parameters.
%          Any fields omitted will be set to default values.
% p_out: A struct containing the derived parameters.

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

implant.IC = 'CIC4';

implant = Ensure_field(implant, 'protocol', 'condensed');
implant = Ensure_field(implant, 'rf_freq',  5e6);	% Hz

% Minimum period in microseconds:
implant.MIN_PERIOD_us = 31.2;

implant.default_phase_width   =  9.6;
implant.default_phase_gap     =  4.8;
