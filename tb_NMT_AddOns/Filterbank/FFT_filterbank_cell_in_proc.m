function out = FFT_filterbank_cell_in_proc(p, cell_in)
% function out = FFT_filterbank_cell_in_proc(p, cell_in)
%
% Quadrature FIR filterbank implemented with FFT.
%
% [out, out_ant] = FFT_filterbank_cell_in_proc(p, cell_in)
%
% p:                    Parameter struct (see below).
% cell_in{1} = audio:   A sampled audio signal.
% cell_in{2,end}:       Outputs of process just before FFT_filterbank_cell_in_proc
%                       In case the previous process is Autoc_opt2_Alt_proc.m, 
%                       cell_in{2,end} corresponds to F0
%
% out{1}:               Frequency-Time-indexed Matrix (FTM) of complex filter outputs.
% out{2}:               F0 [Hz] (only if previous process was Autoc_opt2_Alt_proc.m.
%
% Fundamental parameters:
%   audio_sample_rate:  The sample rate for the audio input signal.
%   analysis_rate:      The number of input blocks analysed per second.
%   window:             The FFT window.
%   block_length:       The number of samples in an input block (FFT length).
% Derived parameters:
%   block_shift:        The number of new samples in each block.
%   bin_freq:           The FFT bin frequency spacing.
%   num_bins:           The number of bins that will be retained.
%   char_freqs:         The centre frequency of each bin.
%
% The analysis_rate is quantised to a sub-multiple of the audio_sample_rate.
%
% There are 3 cases for window and block_length:
% If window is supplied, then block_length is taken from it.
% If block_length is supplied, then a Hann window is used.
% If neither is supplied, then a 128-point Hann Window is used.
%
% Adapted from: FFT_filterbank_proc.m (NMT)
% Adapted by Matthias Milczynski, ExpORL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%    $Revision: 20 $
%        $Date: 2/09/05 4:49p $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0	% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	out = feval(mfilename, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1	% Parameter calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fundamental parameters:

	p = Ensure_field(p, 'audio_sample_rate', 16000);	% Hz
	p = Ensure_field(p, 'analysis_rate',       500);	% Hz

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Optional parameters:

	if isfield(p,'window')
		p.block_length = length(p.window);
	end;

	p = Ensure_field(p,'block_length', 128);
	p = Ensure_field(p,'window',       Cos_window(p.block_length, 'Hann'));
	
	% Set to empty matrix to zero pad at start:
	p = Ensure_field(p, 'buffer_opt', 'nodelay');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Derived parameters:

	p.window_length = length(p.window);
	p.block_shift	= ceil(p.audio_sample_rate / p.analysis_rate);
	p.analysis_rate	= p.audio_sample_rate / p.block_shift;

	p.num_bins		= p.block_length/2 + 1;
	p.bin_freq		= p.audio_sample_rate/p.block_length;
	p.bin_freqs		= p.bin_freq * (0:p.num_bins-1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Display parameters:
	% These fields are not used by this file,
	% and are only set up for GUI_FTM.
	% We are careful not to overwrite any fields that have been set-up
	% by a higher-level function which calls this function.

	p = Ensure_field(p,'char_freqs', p.bin_freqs);
	p = Ensure_field(p,'sample_rate',p.analysis_rate);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	out = p;	% Return parameters.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2	% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    audio = cell_in{1};
	u = buffer(audio, p.block_length, p.block_length - p.block_shift, p.buffer_opt);
	v = u .* repmat(p.window, 1, size(u,2));	% Apply window
	u = fft(v);									% Perform FFT to give Frequency-Time Matrix
	u(p.num_bins+1:end,:) = [];					% Discard the symmetric bins.
    out = {u, cell_in{2}};
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;
