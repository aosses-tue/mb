function uA = Current_level_to_microamps_proc(p, cl)

% Current_level_to_microamps_proc: Convert current level to microamps.
% NaN values represent absent pulses, i.e. current = 0 uA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0  % Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    uA = feval(mfilename, []);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1  % Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = Ensure_implant_params(p);
    p.implant.CURRENT_EXP = log(p.implant.CURRENT_uA_MAX/p.implant.CURRENT_uA_MIN)/p.implant.CURRENT_LEVEL_MAX;
    
    uA = p;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2  % Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check that current levels are within bounds.
    % Note that NaNs do not cause an error.

    if any(cl > p.implant.CURRENT_LEVEL_MAX)
        error('Current level too high.');
    elseif any(cl < 0)
        error('Current level negative.');
    end
    
    uA = p.implant.CURRENT_uA_MIN * exp(p.implant.CURRENT_EXP * cl);

    % Special case for absent pulses:
    empty = isnan(cl);
    uA(empty) = 0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
