function out = FreedomMicResp_proc(p, in)
% function out = FreedomMicResp_proc(p, in)
%
% % Standalone example:
%   p.audio_sample_rate = 16000;
%   p = FreedomMicResp_proc(p);
%
% Programmed by Matthias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case 1
        switch computer
            case 'PCWIN' % Matthias' dir
                root = 'C:\Documents and Settings\u0051992\My Documents\Internal\data';
            case 'GLNX86' % Matthias' dir
                root = fullfile(getenv('INTPATH'), 'data', 'speech', 'noise');
            case 'GLNXA64' % ORL-WRK-0089
                root = '/home/alejandro/';
        end
        dir_old = cd;
        dir_current = get_current_directory(mfilename('fullpath'));
        
        p       = Ensure_field(p, 'noise_file', '../Sounds/noiseFreedom.wav'); % older name: noise1.wav
        
        % 128 filter taps
        p       = Ensure_field(p, 'fl', 128);
        
        try
            cd(dir_current)
            [noise, fs] = wavread(p.noise_file);
            cd(dir_old)
            
            noise   = resample(noise, p.audio_sample_rate, fs);
            p.noise_fs = fs;
            NFFT    = 8192;
            NFFTh   = NFFT/2;
            % chop off 1 sencond at end and beginning
            noise   = noise(fs:end-fs);
            fftresp = fft(noise, NFFT);
            frequencies = 0:1/(NFFTh - 1):1;
            p.b     = firls(p.fl, frequencies, abs(fftresp(1:NFFTh)));
        catch
            error([mfilename '.m: noiseFreedom.wav not found, loading default parameters'])
        end
        out     = p;
    case 2
        out = filter(p.b, 1, in);
end

end