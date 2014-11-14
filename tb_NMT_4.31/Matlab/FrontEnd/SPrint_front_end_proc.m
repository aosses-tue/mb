function y = SPrint_front_end_proc(p, x)

% SPrint_front_end_proc: SPrint analog front end processing (without AGC).
%
% y = SPrint_front_end_proc(p, x)
%
% p: Parameter struct.
% x: Audio input.
% y: Audio output.

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

	y = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Default parameter values:
	
	% Input signal source:
	p = Ensure_field(p, 'source', 'HS8');
	
	% Manual sensitivity setting on SPrint LCD:
	p = Ensure_field(p, 'sensitivity',		12);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% These parameters should almost never be changed:
	
	% From the SPrint hardware specification:
	% +750 mV gives positive full scale (12 bit data)
	p = Ensure_field(p, 'adc_conversion_scaling', 2048/750);
	
	% This value is determined by the SPrint software.
	% It is the same for all strategies.
	p = Ensure_field(p, 'adc_output_gain',	4);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate derived parameters:

	switch p.source
	
	case {'HS8', 'Headset', 'Freedom', 'freedom'}
		p.pre_emphasis      	=  0;
		p.pre_amp_gain_dB 		= 20;
	
    case {'CP810','SP15'}
		p.pre_emphasis      	=  1;
        if ~isfield(p,'pre_amp_gain_dB')
            p.pre_amp_gain_dB 	= -36;
        end
        %%% Added by Alejandro Osses 07.May.13
        
	case {'Lapel microphone', 'Stereo plug ring load'}
		p.pre_emphasis      	=  0;
		p.pre_amp_gain_dB 		= 20;
		
	case {'Personal audio cable', 'TV cable', 'Stereo plug ring open'}
		p.pre_emphasis      	=  0;
		p.pre_amp_gain_dB 		=  0;
		
	case {'Telephone adaptor', 'Stereo plug ring ground', 'Mono plug'}
		p.pre_emphasis      	=  1;
		p.pre_amp_gain_dB 		= 20;	
		
	otherwise	
		error(sprintf('Unknown source: %s', p.source));	
    end
    
    p.sensitivity_gain_dB   = SPrint_sensitivity_gain_dB(p.sensitivity); % > 0
    p.pre_amp_gain          = 10 ^ (p.pre_amp_gain_dB / 20); % > 0
    p = Ensure_field(p, 'reference_mV', 4.7214 / p.pre_amp_gain);
    p.front_end_gain_dB		= p.pre_amp_gain_dB + p.sensitivity_gain_dB;
    p.front_end_gain		= 10 ^ (p.front_end_gain_dB / 20);
	p.front_end_scaling		= p.reference_mV * p.front_end_gain * p.adc_conversion_scaling * p.adc_output_gain;
   
	y = p;	% Return parameters.	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	y = x * p.front_end_scaling;
    
    if p.DEBUG == 1
        
        try
            CP_name = 'CP0101';
            Source = 'NM';
            Y_Limits = [-40 20];
            CP_name = [CP_name '-' p.map_name];
            Controlpoint2Fig(y, CP_name, Source, Y_Limits)
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
