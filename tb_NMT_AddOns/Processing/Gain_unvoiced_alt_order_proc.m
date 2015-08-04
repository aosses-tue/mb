function out = Gain_unvoiced_alt_order_proc(p, cell_in)
% function out = Gain_unvoiced_alt_order_proc(p, cell_in)
%
% Apply gain in dB (multiply by a constant).
% 	v = Gain_proc(p, u)
%
% Adapted from Gain_proc.m, NMT
%
% Programmed by Matthias Milczynski. Comments by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	out = feval(mfilename, []);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p = Ensure_field(p, 'gains_unvoiced_dB', 0.0);
    out = p;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(isfield(p,'DEBUG'))
        fprintf( 1,'Inside: %s\n', mfilename ); 
    end 
    v = cell_in{1,1};
    noAM_idx = cell_in{1,4};
    gains = From_dB(p.gains_unvoiced_dB);
    
	if length(gains) == 1     
        v(:, noAM_idx) = v(:, noAM_idx) * gains;
	else
		v(:, noAM_idx) = v(:, noAM_idx) .* repmat(gains, 1, size(v(:, noAM_idx), 2));
    end
    
	out = {v, cell_in{1,2}, cell_in{1,3}, cell_in{1,4}, cell_in{1, 5}};
    %out = v;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
