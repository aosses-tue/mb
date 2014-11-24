function p = Ensure_implant_params(p)

% Ensure_implant_params: Ensure implant parameters are valid.
%
% p_out = Ensure_implant_params(p_in)
%
% p_in:  A struct containing the clinical parameters.
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
    p = [];
end

p         = Ensure_field(p,         'implant', []);
p.implant = Ensure_field(p.implant, 'IC',      'CIC4'); % Changed to CIC4 from CIC3
p.implant = feval(['Ensure_', p.implant.IC, '_params'], p.implant);

p         = Ensure_field(p, 'phase_width', p.implant.default_phase_width);
p         = Ensure_field(p, 'phase_gap',   p.implant.default_phase_gap);
