function out = CBPF_filterbank_proc(p, in)
% function out = CBPF_filterbank_proc(p, in)
%
% Complex bandpass-filter filterbank formulation of FFT_filterbank proc. 
% Uses linear phase terms (lpt) to frequency-shift the lowpass filter 
% defined by the window function used (usually 128 bin hann window). 
% Works as other NMT procs. 
%
% NOTE: Does not work with windows provided by NMT!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case 0
        out = feval(mfilename, []);
        
    case 1
        p = Ensure_field(p, 'audio_sample_rate', 16000);
        p = Ensure_field(p, 'analysis_rate', 1800);
        p = Ensure_field(p, 'num_bands', 22);
        p = Ensure_field(p, 'downsample', 1);
        p = Ensure_field(p, 'scaleFilters', 0);
        p = Ensure_field(p, 'cbpf_type', 'complex');
        if isfield(p, 'window')
            p.block_length = length(p.window);
        end
        p = Ensure_field(p, 'window_length', 128); 
        p = Ensure_field(p, 'block_length', 128);
        p = Ensure_field(p, 'window', Cos_window(p.window_length, 'Hann'));
        %p = Ensure_field(p, 'window', hanning(p.block_length, 'periodic'));
        p = Ensure_field(p, 'window_type', 'periodic');
        p.block_shift	= ceil(p.audio_sample_rate/p.analysis_rate);
        p.D = p.block_shift;
        p.analysis_rate	= p.audio_sample_rate/p.block_shift;
        p.num_bins		= p.block_length/2 + 1;
        p.bin_freq		= p.audio_sample_rate/p.block_length;
        p.bin_freqs		= p.bin_freq * (0:p.num_bins-1);
        p = Ensure_field(p,'char_freqs', p.bin_freqs);
        p = Ensure_field(p,'sample_rate', p.analysis_rate);
        p = Vector_sum_proc(p);
        h = zeros(p.window_length, p.num_bins);
        h(:, 1) = p.window;
        sh = 1;
        for i = 1:p.num_bins
           lpth = createLPT(sh, i, p.window_length, p.block_length, ...
               p.cbpf_type).'; 
           h(:, i) = lpth.*h(:, 1);
        end
        % filters as used in a CI (based on vector-sum)
        p.hci = p.weights*(h.'); 
        if strcmp(p.cbpf_type, 'complex')
            % Otherwise we get a gain of 6dB (x 2)
            p.hci = p.hci/2;
        end
        out = p;
        
    case 2
        % store impulse responses in columns
        out = cell(2, 1);
        out{1, 1} = in;
        u = zeros(p.num_bands, length(in));
        for i = 1:p.num_bands
            u(i, :) = filter(p.hci(i, :), 1, in); 
        end
        % chop off transient
        chop = p.window_length;
        switch p.window_type
            case 'periodic'
                chop = chop + 1;
        end
        u = u(:, chop:end); 
        % downsample
        if p.downsample
            u = downsample(u.', p.D).';
        end
        %out{2, 1} = u;
        out = u;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end