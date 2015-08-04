function [env,B1,A1,B2,A2] = Extract_MEM_envelope_proc( p, wav )
% function [env,B1,A1,B2,A2] = Extract_MEM_envelope_proc( p, wav )
%
% Extract_MEM_envelope_proc: Extracts envelope of broadband signal for AM
% of the FTM channels.
%
% Reference: Vandali, 2005, JASA
% Programmed by Matthias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 0	% Default parameters
	env = feval(mfilename, []);
    % INIT %
    case 1
        fs = 16e3;
        p = Ensure_field( p, 'fs', fs );
        p = Ensure_field( p, 'lpford', 4);
        p = Ensure_field( p, 'hpford', 2);
        % LPF params
        [lpfB, lpfA] = butter( p.lpford, 300 / ( fs/2 ), 'low' );
        p = Ensure_field( p, 'lpfB', lpfB );
        p = Ensure_field( p, 'lpfA', lpfA );
        % LPF params
        [hpfB, hpfA] = butter( p.hpford, 80 / ( fs/2 ), 'high' );
        p = Ensure_field( p, 'hpfB', hpfB );
        p = Ensure_field( p, 'hpfA', hpfA );
        p = Ensure_field( p, 'lpfSC', 0.5 );
        p = Ensure_field( p, 'bufN', round( 0.012*fs ) );
        env = p;
    % PROCESSING %    
    case 2
        % FWR + LPF
        lpf = filter( p.lpfB, p.lpfA, abs(wav) );
        g  = grpdelay( p.lpfB, p.lpfA );
        L   = length( g );
        % Get allpass filter coeffs for delaying FTM channels that will be
        % modulated to compensate for grpdelay introduced here
        [B1,A1] = iirgrpdelay( p.lpford, 0:1/L:1- 1/L, [0 (1 - 1/L)], g );
        % HPF
        hpf = filter( p.hpfB, p.hpfA, lpf );
        g  = grpdelay( p.hpfB, p.hpfA );
        L   = length( g );
        % Get allpass filter coeffs for delaying FTM channels that will be
        % modulated to compensate for grpdelay introduced here
        [B2,A2] = iirgrpdelay( p.hpford, 0:1/L:1-1/L, [ 0 (1 - 1/L) ], g );
        % HPF
        % Scaling and adding lpf to hpf
        sclpf = lpf * 0.5;
        schpf = hpf + sclpf;
        % HWR
        hwr = 0.5*( abs(schpf) + schpf );
        b = buffer( hwr, p.bufN, 0, 'nodelay' );
        % Normalization
        for j = 1:size( b, 2 )
           b(:,j) = b(:,j)/max(b(:,j)); 
        end
        env = b(1:length(wav));
end
