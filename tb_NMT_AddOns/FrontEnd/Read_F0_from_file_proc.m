function out = Read_F0_from_file_proc(p, wav)
% function out = Read_F0_from_file_proc(p, wav)
%
% Reads pre-calculated F0 values from file. Takes original wav file 
% (e.g. speech in noise) and puts it into output together with read F0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin 
    case 1
        p = Ensure_field(p, 'max_f0', 500);
        p = Ensure_field(p, 'f0var', p.max_f0 - 0.1*p.max_f0 );
        p = Ensure_field(p, 'error_check', 0);
        p = Ensure_field(p, 'removeF0Peaks', 0);
        out = p;
    case 2
        if(p.DEBUG)
            fprintf( 1,'Inside: %s\n', mfilename ); 
        end  
        fid = fopen( p.f0File );
        F0 = [];
        while ~feof(fid)
           val = fgetl( fid );
           val = str2double( val );
		   if p.error_check == 1
			   if val >= p.f0var 
				   % Remove disturbing outliers
				   F0 = [ F0 -1];
			   else
				   F0 = [ F0 val ];
			   end
		   else
			   F0 = [ F0 val ];                   
		   end
        end
        fclose(fid);
        L = length(F0);
        % Remove single local peaks -1 P -1
        if p.removeF0Peaks
            for i = 2:L - 1 
                if F0(i - 1) <= 0 && F0(i + 1) <= 0 && F0(i) > 0
                    F0(i) = -1.0;
                    if p.DEBUG
                        disp('F0 single peak removed');
                    end
                end
            end
        end
        out = { wav, F0 };
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end