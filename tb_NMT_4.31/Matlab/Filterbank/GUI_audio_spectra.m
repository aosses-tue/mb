function GUI_audio_spectra(audio, audio_sample_rate)

% GUI_audio_spectra: GUI to display and play an array of audio spectra.
%
% GUI_audio_spectra(audio, audio_sample_rate)
%
% audio:                An array of audio signals.
%                       Audio can run across rows or columns.
% audio_sample_rate:	The sample rate of each signal.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate spectra:

[num_rows, num_cols] = size(audio);
if (num_rows > num_cols)
	u.audio = audio';
else
	u.audio = audio;
end
[num_audio_signals, len] = size(u.audio);
if (num_audio_signals == 1)
	error('Must have more than one band');
end

u.audio_sample_rate = audio_sample_rate;	

p = FFT_filterbank_proc;
for n = 1:num_audio_signals
	cmplx = FFT_filterbank_proc(p, u.audio(n,:));
	u.mag{n} = sqrt(abs(cmplx));
end

GUI_audio_FTM(u);

