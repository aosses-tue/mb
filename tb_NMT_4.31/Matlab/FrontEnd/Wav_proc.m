function x = Wav_proc(p, a)

% Wav_proc: Read wav file and resample if necessary.
% x = Wav_proc(p, a)
% If the arg is a string, assume it is the name of a .wav file,
% and read it from disk. If not, just pass it through.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent start_index
persistent end_index
persistent num_samples

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	x = feval(mfilename, []);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = Ensure_field(p, 'audio_sample_rate', 16000);
	p = Ensure_field(p, 'wav_sample_rate_tolerance', 1.05);
	x = p;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~ischar(a)
		x = a;		% pass through
		
	else	% Assume it is the name of a wave file
	
		if ~isfield(p, 'chunk_size')
		
			% Read entire file as a single chunk:
			[x, fs] = wavread(a);
			% Allow small sample rate differences:
			rate_ratio = fs / p.audio_sample_rate;
			if (rate_ratio   > p.wav_sample_rate_tolerance)...
			|  (1/rate_ratio > p.wav_sample_rate_tolerance)
				x = resample(x, p.audio_sample_rate, fs);
			end

		else
			
			if isempty(start_index)
				% First chunk:
				size = wavread(a, 'size');
				num_samples = size(1);
				start_index = 1;
				end_index   = p.chunk_size;
			else
				start_index = start_index + p.chunk_size;
				end_index   = end_index   + p.chunk_size;
			end

			if end_index > num_samples
				end_index = num_samples;
			end

			if start_index <= num_samples
				% Read next chunk from file:
				x = wavread(a, [start_index, end_index]);
			else
				% Finished reading file:
				x = [];
				start_index = [];
				end_index   = [];
				num_samples = [];
			end
		end
    end
    
    if isfield(p,'calibrate_microphone')
        if p.calibrate_microphone == 1
            display(['Front end process temporarily added to ' mfilename]);
            ch_b = [1.407226562 4.21434833817899 4.21434833817899 1.407226562];
            ch_a = [2 4.61172342112964 3.66191296364167 0.96951341558667];
            calTaps_Fmic = Calibrate_xPC4_Microphone_Response;
            x = filter(ch_b, ch_a, x);
            x = filter(calTaps_Fmic, 1, x);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
