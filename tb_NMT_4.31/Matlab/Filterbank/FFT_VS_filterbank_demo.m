function FFT_VS_filterbank_demo(audio)

% FFT_VS_filterbank_demo: GUI demonstrating FFT VS filterbank.
% See GUI_audio_spectra, GUI_audio_FTM for more information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(audio)
	audio = wavread(audio);
end

p1.num_bands = 6;
p1.audio_sample_rate = 16000;
p1.analysis_rate = p1.audio_sample_rate;
p1 = Append_process(p1, 'FFT_filterbank_proc');
p1 = Append_process(p1, 'Vector_sum_proc');

v1 = Process(p1, audio);
GUI_FTM(p1, v1, 'Filterbank output');

audio = [audio'; real(v1)];

GUI_audio_spectra(audio, p1.audio_sample_rate);
