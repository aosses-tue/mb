function [name path extension] = Split_file_name(filename_w_path)
% function [name path extension] = Split_file_name(filename_w_path)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       filename_w_path = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-calibrated\modus-1_v2-2filt.wav';
%       [name path extension] = Split_file_name(filename_w_path);
%       disp(['file name     : ' name])
%       disp(['file extension: ' extension])
%       disp(['Located at    : ' path])
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 1/7/2014
% Last update: 1/7/2014 % Update this date manually
% Last used: 1/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[path,name,extension] = fileparts(filename_w_path);

path = [path delim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end