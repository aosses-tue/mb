function [v, sub, sat] = LGF_proc(p, u)

% LGF_proc: Loudness Growth Function
% function [v, sub, sat] = LGF_proc(p, u)
%
% Inputs:
% p:            Parameter struct, with the following fields:
% p.Q:            Percentage decrease of output when input is 10 dB below sat_level.
% p.lgf_alpha:    Curve shape factor.
% p.base_level:   Input value which maps to 0.
% p.sat_level:    Input value which maps to 1.
% p.sub_mag:      Output value used for inputs less than base_level (negative or zero).
% u:            Input magnitude vector or FTM.
%
% Outputs:
% v:            Magnitude in range 0:1 (proportion of dynamic range).
% sub:          Logical FTM indicating the values that were below base_level.
% sat:          Logical FTM indicating the values that were above sat_level.

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

	v = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Defaults:
    p = Ensure_field(p,'scale_LGF_input',0); % Added by Alejandro Osses
    
	p = Ensure_field(p,'base_level',  4/256);
	p = Ensure_field(p,'sat_level', 150/256);
	p = Ensure_field(p,'Q',          20);
	p = Ensure_field(p,'sub_mag',    -1e-10);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Derived parameters:
	
	if (p.base_level > 0)
		% for information, not used in processing.
		p.lgf_dynamic_range = 20*log10(p.sat_level/p.base_level); 
    end
    
    if ~isfield(p,'lgf_alpha')
        p.lgf_alpha	= LGF_alpha(p.Q, p.base_level, p.sat_level);
    end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	v = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isfield(p,'scale_LGF_input')
        if p.scale_LGF_input == 1
            u = p.gain_LGF_input * u;
        end
    end
    
	% Scale the input between base_level and sat_level:
	r = (u - p.base_level)/(p.sat_level - p.base_level);

	% Find all the inputs that are above sat_level (i.e. r > 1) 
	% and move them down to sat_level:
	sat = r > 1;		% This is a logical matrix, same size as r. 
	r(sat) = 1;

	% Find all the inputs that are below base_level (i.e. r < 0) 
	% and temporarily move them up to base_level:
	sub = r < 0;		% This is a logical matrix, same size as r.
	r(sub) = 0;

	% Logarithmic compression:
	v = log(1 + p.lgf_alpha * r) / log(1 + p.lgf_alpha);

	% Handle values that were below base_level:
	v(sub) = p.sub_mag;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

