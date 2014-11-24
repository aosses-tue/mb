function implant = Ensure_CIC1_params(implant)

% Ensure_CIC1_params: Ensure implant parameters are valid for CIC1.
%
% p_out = Ensure_CIC1_params(p_in)
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

implant.IC = 'CIC1';

implant = Ensure_field(implant, 'protocol', 'expanded');
implant = Ensure_field(implant, 'rf_freq',  2.5e6);   % Hz

implant.CURRENT_LEVEL_MAX	= 239;

implant.CURRENT_uA_MIN		=   19.63;
implant.CURRENT_uA_MAX		= 1599.00;
