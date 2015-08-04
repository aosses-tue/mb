function out = Cos_ramp_proc(p, in)
% function out = Cos_ramp_proc(p, in)
%
% Apply on and off ramp to input.
%
% Programmed by Matthias Milczynski, comments by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        out = feval(mfilename, []);
        
    case 1
        
        p = Ensure_field(p, 'ramp_sample_rate', p.audio_sample_rate);
        % attack-time in ms
        p = Ensure_field(p, 'ramp_atime', 15);
        % release-time in ms
        p = Ensure_field(p, 'ramp_rtime', p.ramp_atime);
        out = p;
    case 2
		if isstruct(in)
			x = sum(in.ftmo)';
		else
			x = in;
		end
        ramp = cos_ramp(p.ramp_sample_rate, length(x), p.ramp_atime, ...
            p.ramp_rtime);
        out = ramp(:).*x(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end