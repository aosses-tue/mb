function out = CP810MicResp_proc( p, in )
% function out = CP810MicResp_proc( p, in )
%
% Emulates processing for CP810 microphones (SP15 Speech Processors)
%
% Champ-filter
% Frequency response compensation
% Directionality filter (without Beamforming)
%
% Programmed by Alejandro Osses, ExpORL, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1	% Parameter calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        mic_model = ['xPC' num2str(p.micfront - 15000)];
        
        % champ response for both channels (Front and Rear) (IIR Filter)
        p.ch_b = [1.407226562 4.21434833817899 4.21434833817899 1.407226562];
        p.ch_a = [2 4.61172342112964 3.66191296364167 0.96951341558667];
        
        % FIR Filter
        [calTaps_Fmic, calTaps_Rmic] = Calibrate_SP15_Microphone_Response(mic_model,p.audio_sample_rate );

        p.b_Front   = calTaps_Fmic; % numerator coefficients
        p.b_Rear    = calTaps_Rmic;
        
        % Equivalent to 'SP15R1 Directionality' Simulink block

        directionality = 2; % from 1 to 4
        disp(['Inside ' mfilename ' Loading default directionality...to fix this in the future.'])
        
        switch directionality
            case 1 % Omni front
                directionality = 'Front Omni';
            case 2 % Near Omnidirectional
                directionality = 'Zoom Disabled';
            case 3 % Zoom
                directionality = 'Zoom On';
            case 4 % Flat front
                directionality = 'Front Flat';
            case 5 % Passthrough front
                directionality = 'Front Passthrough';
        end

        shape_select = directionality;
        beamTaps = load('beamformerFilters.mat');
        ssel = shape_select;
        ssel(ssel == ' ') = '_';

        fmic_taps = beamTaps.(ssel).front;
        rmic_taps = beamTaps.(ssel).rear;

        p.fmic_taps = fmic_taps;
        p.rmic_taps = rmic_taps;
        
        out = p; % Modified structure
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2	% Processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(p.DEBUG)
            fprintf( 1,'Inside: %s\n', mfilename ); 
        end  
        % Equivalent to 'CP810 Calibration Filter - ParMate' Simulink block
        inChamped   = filter(p.ch_b     , p.ch_a, in);
        out2Directionality = filter(p.b_Front  , 1, inChamped);
        
        % Equivalent to 'SP15R1 Directionality' Simulink block
        out = filter(p.fmic_taps, 1, out2Directionality);
        
        if p.DEBUG == 1
            try
                CP_name = 'CP0201';
                Source = 'NM';
                Ymin = 0;
                Y_Limits = [Ymin Ymin + 60];
                CP_name = [CP_name '-' p.map_name];
                Controlpoint2Fig(out, CP_name, Source, Y_Limits)
            end
        end
    end
end
