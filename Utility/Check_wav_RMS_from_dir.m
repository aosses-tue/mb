function Check_wav_RMS_from_dir(dir)
% function Check_wav_RMS_from_dir(dir)
%
% 1. Description:
%       Gives RMS values [dBFS] of all wav files inside 'dir'
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       Check_wav_RMS_from_dir; % then select folder containing wav files
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 20/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    dir = uigetdir(Get_TUe_paths('outputs'));
    dir = [dir delim];
end

filenames = Get_filenames(dir,'*.wav');

if length(filenames) ~= 0
    
    disp([mfilename '.m: RMS values of file inside ' dir])
    
    for i = 1:length(filenames)
        x = Wavread([dir filenames{i}]);
        fprintf('rms value = %.3f [dBFS], filename: %s\n',rmsdb(x),filenames{i})
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
