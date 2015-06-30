function out = AM_BS_corr_proc(p, cell_in)

% function out = AM_BS_corr_proc(p, cell_in)
%
% Amplitude modulation of filterbank channels with extracted F0
% values. 
%
% Inputs:
% map:      Parameter struct
% input:    Cell array two cells (Processing only):
%                           - FTM of filterbank channels (1)
%                           - F0 values of current audio file (2)
%                           - AM_idx (3)
%                           - noAM_idx (4)
%
% Outputs:
% out:     Parameter struct (Initialization) or amplituce modulated FTM 
%          (Processing) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: K.U.Leuven, Lab. Exp. ORL, Matthias Milczynski
%        $Date: 04/09/2006 $
%      Authors: Matthias Milczynski
% Edited by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        return;
    case 1
    %%%%%%%%%%%%%%%%
    % Initialization
    %%%%%%%%%%%%%%%%
        
        p = Ensure_field(p, 'audio_sample_rate' , 16000);
        p = Ensure_field(p, 'nam_bs_outs'       , 1);
        p = Ensure_field(p, 'bApplyModulator'   , 1);
        out = p;
        
    case 2
    %%%%%%%%%%%%
    % Processing
    %%%%%%%%%%%%
        if(isfield(p,'DEBUG') ) 
            fprintf( 1,'Inside: %s\n', mfilename ); 
        end 
        
        penv        = cell_in{1,1}; % FTM of filterbank channels
        Modulator   = cell_in{1,2}; % Modulator values of current audio file
        AM_idx      = cell_in{1,3}; % AM_idx   - Amplitude Modulation Indexes
                                % noAM_idx - Indexes of samples with no AM (not used in this procedure) 
        rms_ACE     = cell_in{1,5};
        tstep       = 1/p.analysis_rate;

        % Get AM segments    
        d       = diff( AM_idx );
        bnds    = find( d > 1 );
        Modulator_New  = Modulator;
        st = 1;
        s = zeros(1, size(penv, 2 ));
        if ~isfield(p,'cisim') % AOV, Originally: if ~p.cisim
            for i=1:length( bnds )
                en = bnds(i);
                idx = AM_idx( st:en );
                t = 0:tstep:length( Modulator( idx ) )/p.analysis_rate - tstep;
                mpeak = max(max(penv(:, idx )));
                Modulator_New(idx) = scsin( mpeak, Modulator( idx ), p.base_level, p.sat_level);
                st = en + 1;
            end

            if( st < size( penv, 2 ) )
                idx = AM_idx( st:end );
                mpeak = max(max(penv(:, idx )));
                t = 0:tstep:length( Modulator( idx ) )/p.analysis_rate - tstep;
                if length(idx) ~= 0;
                    Modulator_New(idx) = scsin( mpeak, Modulator( idx ),p.base_level, p.sat_level);
                end
            end

            if(p.plot_modulator)
                figure;
                t = 0:tstep:length(Modulator)/p.analysis_rate - tstep;
                plot(t, Modulator_New, t, Modulator);grid on; legend('Modulator New','Modulator')
            end
        end
        for i = 1:size(penv, 2)
            if p.bApplyModulator == 1
                penv(:, i) = penv(:, i).*Modulator_New(i);
            end
        end
        if isfield(p,'rmsequalize') 
            for j=1:p.num_bands
               penv(j, :) = eq_rms(penv(j, :), rms_ACE(j)); 
            end
        end
        if p.nam_bs_outs == 1
            out = penv;
        else
            out = {penv, cell_in{1,6}};
        end
end     
