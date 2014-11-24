function seq = Read_sequence(file_name)
% Read_sequence: Reads a stimulus sequence file.
% function seq = Read_sequence(file_name)
% If an output arg is omitted, uses the the file base name as the workspace name.
% The following file extensions are supported:
%
%    Ext    Type    Contents    Tool
%   ----    ------  --------    ----------------------------------
%   .rft    text    frames      Output of RxFrames.exe
%   .rfd    text    frames      Obsolete
%   .quf    binary  frames      Input to SPrint streaming software.
%
% Examples:
%   s1 = Read_sequence('foo.rft');  % creates s1
%   s2 = Read_sequence('foo.rfd');  % creates s2
%   s3 = Read_sequence('foo.quf');  % creates s3
%   s4 = Read_sequence('foo');      % reads foo.quf, creates s4
%   Read_sequence('foo.quf');       % creates foo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 1)
    help(mfilename);
    error('File name must be supplied');
end

[path_name, base_name, extension] = fileparts(file_name);

if (isempty(extension))
    % Use default extension:
    extension = '.quf';
    file_name = fullfile(path_name, [base_name, extension]);    
end

switch extension

    case '.rfd'
        seq = Read_RFD(file_name);
        
    case '.rft'
        seq = Read_RFT(file_name);

    case '.quf'
        seq = Read_QUF(file_name);
        
    otherwise
        error('Unknown file extension');
end

if (nargout == 0)
    assignin('caller', base_name, seq); % assign in caller's scope.
    seq = [];                           % don't need it in ans too.
end
