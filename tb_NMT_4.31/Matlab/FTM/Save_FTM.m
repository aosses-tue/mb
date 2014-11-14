function Save_FTM(ftm, file_name)
% Save_FTM: Saves a Frequency-Time Matrix to a binary file.
% function function Save_FTM(ftm, file_name)
% If file_name is omitted, uses the ftm argument's workspace name as the file base name.
% Examples:
%	Save_FTM(out, 'foo.mat');	% creates foo.mat
%	Save_FTM(out, 'foo');		% creates foo.mat
%	Save_FTM(out);				% creates out.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(1,2,nargin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine file name and type:

path		= [];
base_name	= [];

if (nargin == 2)
	% file_name argument was supplied.
	[path, base_name, extension] = fileparts(file_name);
	if (~isempty(extension) & ~isequal(extension, '.mat'))
		error('Only .mat is supported');
	end
end

if (isempty(base_name))
    base_name = inputname(1);
	if (isempty(base_name))
		error('Must supply a file name');
	end
end

full_name = [path, base_name];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale data to uint16 and save:

max_value = 65535;

m = max(max(ftm));
if (round(m * max_value) > max_value)
	warning(sprintf('max value %g exceeds 1', m));
end

ftm16 = uint16(ftm * max_value + 0.5);

save(full_name, 'ftm16');