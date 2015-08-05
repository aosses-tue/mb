function out = ExpORLResynth_proc(p, in)
% function out = ExpORLResynth_proc(p, in)
%
% Vocoder based CI simulation resynthesizer proc written in NMT style.
% Approaches were mainly developed in Laneau 2005, JASA.
%
% Crucial parameters:
%
% p.resynth_carrier : 'sinus', 'noise_preFilt' or 'noise_postFilt'
% p.resynth_fdes : 'IIR' or 'FIR'
% p.resynth_fbank : 'CISIM' (Laneau, 2005) or 'ACE'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        out = feval(mfilename, []);
        
    case 1
        % TODO: Zero-padding before filtering to avoid chopped off ending
        p = Ensure_field(p, 'audio_sample_rate', 16e3);
        p = Ensure_field(p, 'resynth_rate', p.audio_sample_rate);
        p = Ensure_field(p, 'resynth_carrier', 'noise');
        p = Ensure_field(p, 'resynth_order', 'postFilt');
        % filter order (required for IIR case, FIR will be determined automatically)
        p = Ensure_field(p, 'resynth_ford', 4);
        % filter type 
        p = Ensure_field(p, 'resynth_ftype', 'IIR');
        % use ACE filterbank or CISIM (greenwood based) filterbank
        p = Ensure_field(p, 'resynth_fbank', 'ACE');  
        switch p.resynth_fbank
            case 'ACE_fitted'
                [p.resynth_B, p.resynth_A] = eval(['Create_' p.resynth_ftype ...
                    '_resynth_filters(p);']);
            case 'ACE'
                cbpf = [];
                cbpf.num_bands = p.num_bands;
                cbpf = Ensure_field(cbpf, 'cbpf_type', 'real');
                cbpf = CBPF_filterbank_proc(cbpf);
                p.resynth_ford = 128;
                p.resynth_B = cbpf.hci;
                p.resynth_A = ones(p.num_bands, 1);
        end
        if strcmp(p.resynth_ftype, 'FIR')
           % we only consider FIR, order is p.resynth_ord + 1, thus
           % grpdelay is resynth_ford / 2; (N - 1)/2
           p.gd = p.resynth_ford/2;
        end
        p = Ensure_field(p, 'sidebands', 2);
        out = p;
        
    case 2
        % Resample FTM
        ind = isnan(in);
        in(ind) = 0;
        [n, d] = rat(p.resynth_rate/p.analysis_rate);
        ftm = resample(in', n, d)'; 
        if strcmp(p.resynth_ftype, 'FIR') && strcmp(p.resynth_order, ...
            'postFilt')
            ftm = [ftm zeros(p.num_bands, p.gd)];        
        end
        switch p.resynth_carrier
            case 'noise'
                carrier = randn(size(ftm));
        end
        switch p.resynth_order    
            case 'postFilt'
                % First apply envelope to carrier then filter. This should
                % be done in F0mod since modulating with a HWR envelope can
                % alter the spectrum of the noise band. 
                ftmo = carrier.*ftm;
                for c = 1:size(ftm, 1)
                   ftmo(c, :) = filter(p.resynth_B(c, :), ...
                       p.resynth_A(c, :), ftmo(c, :)); 
                end
                if strcmp(p.resynth_ftype, 'FIR')
                    ftmo = ftmo(:, p.gd+1:end);
                end
            case 'preFilt'
                % first filter carrier and then modulate with envelopes
                % should be used for F0mod with FIR filter (?). 
                for c = 1:size(ftm, 1)
                   carrier(c, :) = filter(p.resynth_B(c, :), ...
                       p.resynth_A(c, :), carrier(c, :)); 
                end
                ftmo = carrier.*ftm;
        end
        %out = sum(ftmo, 1);
		out = [];
        out.ftmo=ftmo;
        out.ftm=ftm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end