function out = LPF_zero_ph_proc(p, cell_in)
% function out = LPF_zero_ph_proc(p, cell_in)
%
% LPF_AM_proc: Lowpass filtering and amplitude modulation of
% power-sum_envelope
%
% INPUTS:
%   p       :   map struct
%   cell_in :   cell_in{1,1} contains FTM, 
%               cell_in{1,2} contains F0 (one F0 per F0-estimator buffer length)
%
% OUTPUTS:
% out:     Parameter struct (Initialization) or LPF and amplitude modulated FTM 
%          (Processing) 
%   out{1}:
%   out{2}:     Modulator
%   out{3}:     AM_idx      - indexes of segments to be amplitude modulated (voiced segments)
%   out{4}:     noAM_idx    - indexes of unvoiced segments
%   out{5}:
%   out{6}:     F0
%
%      Matthias Milczynski 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0	% Default parameters
        out = feval(mfilename, []);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1 % Initialization
        p = Ensure_field(p, 'audio_sample_rate', 16000);
        p = Ensure_field(p, 'bApplyLPF_IIR', 1); % Added by AO
        out = p;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % Processing
        
         if(isfield(p,'DEBUG')) 
              fprintf( 1,'Inside: %s\n', mfilename ); 
          end 

        % Get length of F0 array and FTM
        v           = cell_in{1,1};
        F0          = cell_in{1,2};
       
        F0_len      = size(F0, 2);
        pEnv_len    = size(v, 2);
        v_prev      = v;
        
        % Time scale unit for F0_in in relation to pEnv 
        ac_ts       = F0_len/pEnv_len;          
        f_update    = p.analysis_rate;

        F0_mat      = repmat(F0(:), 1, round(1/ac_ts));
        F0_mat      = F0_mat';
        F0_in       = F0_mat(:); % F0_in would be uniquely a column vector

        if( size(v, 2) > size(F0_in, 1))    % Size adjustment, it is needed:
                                            %   Num of raws in F0_in =  Num of columns in v
            delta   = size( v, 2 ) - size( F0_in, 1 );
            F0_in   = [F0_in;ones(delta, 1)];
            
        elseif( size( v, 2 ) < size( F0_in, 1 ) )
            
            delta   = size( F0_in, 1) - size( v, 2 );
            F0_in   = F0_in( 1:end-delta, 1 );    
            
        end
        
        if(size(v, 2) ~= size(F0_in, 1))
            error('!F0 vector and pEnv row-length still not equal!');
        end
        
        % Amplitude modulation
        noAM_idx    = find(F0_in <= 0);
        AM_idx      = find(F0_in > 0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAKE SURE EACH MODULATION SEQUENCE STARTS AT PHASE 0
        
        % collect connected voiced regions
        d           = diff( AM_idx );
        bnds        = find(d > 1);
        
        st          = 1;
        s1          = zeros( 1, size( v, 2 ) );
        s2          = zeros( 1, size( v, 2 ) );
        
        for i =1:length( bnds )
            en      = bnds(i);    
            switch p.modshape
                case 'sine'
                    s1( 1, AM_idx( st:en ) ) = sin(      2*pi / f_update * cumsum(   F0_in( AM_idx( st:en ) ) )   )';
                case 'sawsharp'
                    s1( 1, AM_idx( st:en ) ) = sawtooth( 2*pi / f_update * cumsum( 2*F0_in( AM_idx( st:en ) ) ),0 )';
                    s2( 1, AM_idx( st:en ) ) = sawtooth( 2*pi / f_update * cumsum(   F0_in( AM_idx( st:en ) ) ),0 )';                    
                otherwise
                    error('Unknown modulation shape p.modshape');
            end        
            st = en + 1;
        end
        
        if( st < size( v, 2 ) )
            switch p.modshape
                case 'sine'
                    s1( 1, AM_idx( st:end ) ) = sin(      2*pi / f_update * cumsum(   F0_in( AM_idx( st:end ) ) )   )';
                case 'sawsharp'
                    s1( 1, AM_idx( st:end ) ) = sawtooth( 2*pi / f_update * cumsum( 2*F0_in( AM_idx( st:end ) ) ),0 )';
                    s2( 1, AM_idx( st:end ) ) = sawtooth( 2*pi / f_update * cumsum(   F0_in( AM_idx( st:end ) ) ),0 )';                    
                otherwise
                    error([mfilename '.m: Unknown modulation shape p.modshape']);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        F0_in(noAM_idx) = 0;
        switch p.modshape
            case 'sine'
                Modulator = 0.5 + 0.5*s1;
            case 'sawsharp'
                s1  = s1*0.5 + 0.5;
                idx = s1 <= p.base_level;
                s1(idx) = p.base_level + 0.001;  
                s2 = 1/2*(s2 + abs(s2));
                idx = s2 > 0;
                s2( idx ) = 1;
                figure;
                plot(1:length(s1),s1,1:length(s2),s2);grid on;
                s1 = s1.*s2;
                Modulator = s1;
             otherwise
                error([mfilename '.m: Unknown modulation shape p.modshape']);  
        end
        Modulator(noAM_idx) = 1;
        
        if length(AM_idx) ~= 0 % Then voiced segments were found
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Only apply LPF to voiced parts
            v_AM            = v(:, AM_idx);        
            [B, A]          = butter(8, 60/(p.sample_rate/2), 'low');
            try
                if p.bApplyLPF_IIR
                    vfilt       = filtfilt(B, A, v_AM')';
                else
                    vfilt       = v_AM;
                end
            catch
                vfilt       = v_AM;
                disp([mfilename '.m: Filtering not done: perhaps there was not detected a real voice region']);
            end
            v(:, AM_idx)    = vfilt;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            maxace          = max(v_prev(:, AM_idx), [], 2);
            maxf0m          = max(     v(:, AM_idx), [], 2);
            chRatios        = maxace(:)./maxf0m(:);
            % for i = 1:size(v(:, AM_idx), 2)
            %   v(:, AM_idx(i)) = v(:, AM_idx(i)).*chRatios;
            % end
        end
        
        rmsace          = rmsdb(v_prev');

        out = {v, Modulator, AM_idx, noAM_idx, rmsace, F0_in};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end