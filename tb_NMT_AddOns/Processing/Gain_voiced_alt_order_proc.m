function out = Gain_voiced_alt_order_proc(p, cell_in)
% function out = Gain_voiced_alt_order_proc(p, cell_in)
%
% Edited from v = Gain_proc(p, u)
%
% cell_in{1}    - v - Array: Num_Channels x Num_Time_Slots
% cell_in{2}
% cell_in{3}    - AM_idx - Amplitude modulation index
% cell_in{4}
% cell_in{5}
% cell_in{6}
% out           - Amplified output ( v*Gain )
%
% Dependencies:
%   From_dB (Nucleus MATLAB TB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original script by:
%    Copyright: Cochlear Ltd
%    $Revision: 4 $
%    $Date    : 19/01/05 3:12p $
%    Authors  : Brett Swanson
%
% Edition for F0mod adaption: Matthias Milczynski
% Revision: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	out = feval(mfilename, []);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p = Ensure_field(p, 'gains_voiced_dB', 0.0);
    out = p;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( isfield(p,'DEBUG') )
        fprintf( 1,'Inside: %s\n', mfilename ); 
    end 
    
    v = cell_in{1,1};
    AM_idx = cell_in{1,3};
    gains = From_dB(p.gains_voiced_dB);
        
	if length(gains) == 1     
        v(:, AM_idx) = v(:, AM_idx) * gains;
	else
		v(:, AM_idx) = v(:, AM_idx) .* repmat(gains, 1, size(v(:, AM_idx), 2));
    end
    
    out = {v, cell_in{1,2}, cell_in{1,3}, cell_in{1,4}, cell_in{1, 5}};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
