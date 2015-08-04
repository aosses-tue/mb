function out = Normalize_wav_proc( p, wav )
% function out = Normalize_wav_proc( p, wav )
%
% Peak-normalizes input wav file
%
% Programmed by Matthias Milczynski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case 1
        
        out = p;
         
    case 2
        if(p.DEBUG)
            fprintf( 1,'Inside: %s\n', mfilename ); 
        end  
        out = wav / max( abs( wav ) );
end
