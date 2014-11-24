function [result,sdif,smax,savg] = Compare_sequences(s1, s2, verbose, tolerance)
% Compare_sequences: Compares two sequences
% function [result,sdif,smax,savg] = Compare_sequences(s1, s2, verbose)
%	s1:			sequence 1.
%	s2:			sequence 2.
%	verbose:	flag controlling output display:
%		= 0		No output (default)
%		= 1		Displays messages, max and average differences
%		= 2		Displays messages, and difference sequence
%   tolerance:  parameter tolerance structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(2,4,nargin));
if (nargin < 3)
    verbose = 0;
end
if (nargin < 4)
    tolerance = [];

    % Data fields should match exactly:
    tolerance.channels = 0;
    tolerance.magnitudes = 0;
    tolerance.electrodes = 0;
    tolerance.modes = 0;
    tolerance.current_levels = 0;

    % Timed fields should be almost equal:
    tolerance.phase_widths = 1;
    tolerance.phase_gaps = 2;
    tolerance.periods = 4;
end

% If arguments are strings, assume they are file names:

if isstr(s1)
	s1 = Read_sequence(s1);
end

if isstr(s2)
	s2 = Read_sequence(s2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames1 = fieldnames(s1);
fnames2 = fieldnames(s2);

required_fields = { 'electrodes', 'modes', 'current_levels', 'phase_widths', 'phase_gaps', 'periods' };

if (~all(isfield(s1, required_fields)) | ~all(isfield(s2, required_fields)))
	result = 0;
	if verbose disp('Different fields'); end
	return;
end

n1 = Get_num_pulses(s1);
n2 = Get_num_pulses(s2);

if (n2 ~= n1)
	result = 0;
	if verbose fprintf('Different number of stim pulses: %d %d\n', n1, n2); end
	return;
end

sdif = [];
smax = [];
savg = [];

for n = 1:length(fnames1)
	name = fnames1{n};
	v1 = getfield(s1, name);
	v2 = getfield(s2, name);
	vdif = v1 - v2;			% OK if one or both have length 1
	sdif = setfield(sdif, name, vdif);
	smax = setfield(smax, name, max(abs(vdif)));
	savg = setfield(savg, name, mean(vdif));
end

if verbose == 2
	Disp_sequence(sdif);	% This has N pulses
end
if verbose
	% These have only one pulse
	smax
	savg
end

result = [];

% Compare each of the field present, using the specified tolerance.
if isfield(smax, 'channels')
    result(end+1) = isequal(smax.channels, tolerance.channels);
end
if isfield(smax, 'magnitudes')
    result(end+1) = isequal(smax.magnitudes, tolerance.magnitudes);
end
if isfield(smax, 'electrodes')
    result(end+1) = isequal(smax.electrodes, tolerance.electrodes);
end
if isfield(smax, 'modes')
    result(end+1) = isequal(smax.modes, tolerance.modes);
end
if isfield(smax, 'current_levels')
    result(end+1) = isequal(smax.current_levels, tolerance.current_levels);
end
if isfield(smax, 'phase_widths')
    result(end+1) = (smax.phase_widths <= tolerance.phase_widths);
end
if isfield(smax, 'phase_gaps')
    result(end+1) = (smax.phase_gaps <= tolerance.phase_gaps);
end
if isfield(smax, 'periods')
    result(end+1) = (smax.periods <= tolerance.periods);
end

result = all(result);
