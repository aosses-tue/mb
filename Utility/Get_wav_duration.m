function dur = Get_wav_duration(filename)
% function dur = Get_wav_duration(filename)
%
% 1. Description:
%       Returns total length of an audio file in seconds
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       paths = Get_TUe_paths('db_fastl2007');
%       filename = [paths 'track_38.wav'];
%       Get_wav_duration(filename);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x fs] = Wavread(filename);

dur = max( (1:length(x))/fs );

if nargout == 0
    fprintf('wav file duration: %.3f [s], RMS = %.4f [dBFS]\n',dur,rmsdb(x))
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
