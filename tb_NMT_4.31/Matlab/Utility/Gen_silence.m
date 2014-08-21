function silence = Gen_silence(duration, sample_rate)
% Gen_silence: Generate a sampled-audio silence.
% function silence = Gen_silence(duration, sample_rate)
% duration:		the duration of the silence, in seconds.
% sample_rate:	the audio sample rate, in Hertz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = round(duration * sample_rate);	% Number of time samples.
silence = zeros(N,1);				% Column vector.
