function ftm = Read_FTM(file_name)
% Read_FTM: Reads a Frequency-Time Matrix file.
% function ftm = Read_FTM(file_name)
% If an output arg is omitted, uses the file base name as the workspace name.
%
% Examples:
%	m1 = Read_FTM('foo.mat');	% creates m1
%	m2 = Read_FTM('foo');		% creates m2
%	Read_FTM('foo');			% creates foo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(1,1,nargin));

[path, base_name, extension] = fileparts(file_name);
if (~isempty(extension) & ~isequal(extension, '.mat'))
	error('Only .mat is supported');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data and scale to [0,1]

max_value = 65535;
s = load(file_name);	% s is a struct

ftm = double(s.ftm16) / max_value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If an output arg is omitted, use the file base name as the variable name.

if (nargout == 0)
    assignin('caller', base_name, ftm);
    ftm = [];	% otherwise ans has a copy of ftm too.
end
