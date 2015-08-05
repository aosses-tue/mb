function out = Scale_magnitude_proc(p, v)
% function out = Scale_magnitude_proc(p, v)
%
% Scale_magnitude_proc: Apply scale factor to all magnitudes of a channel
% magnitude sequence
%
% Programmed by Matthias Milczynski. Comments by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
	if isstruct(p) == 0
		return;
	end
end

switch nargin

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0	% Default parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        q = feval(mfilename, []);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1	% Parameter calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         p   = Ensure_field( p, 'magnitude_scaling', 1);
         out = p;    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2	% Processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if( ~isfield( v, 'magnitudes' ) )
            error( 'No magnitudes field in channel magnitude structure' );
        end

        if length( p.magnitude_scaling ) == 1     
            v.magnitudes = v.magnitudes * p.magnitude_scaling; 
        else
            for i = 1:length( v.channels )
                if( v.channels( i ) > 0 )
                    v.magnitudes( i ) = v.magnitudes( i )*p.magnitude_scaling( v.channels( i ) );
                end
            end
        end
        
        out = v;
end
    
    
