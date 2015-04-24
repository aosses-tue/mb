%--------------------------------------------------------------------------
% CaspPreProcInit.m
%--------------------------------------------------------------------------
% calculates filter coefficients and stuff for the basilar-membrane and
% modulation filterbank
%
% usage
%   [BM MF Lp] = CaspPreProcInit(BM, MF, fs, fsIR)
%
% input
%   BM          : structure containing filter coefficients and stuff for
%                 the basilar-membrane filterbank
%   MF          : structure containing filter coefficients and stuff for
%                 the modulation filterbank
%   fs          : sample rate input signal
%   fsIR        : downsampled sample rate of the input signal
%
% output
%   BM          : see above, plus center frequencies are added
%   MF          : see above, plus center frequencies are added
%   Lp          : structure containing filter coefficients for 1 kHz
%                 lowpass for the hair cell transformation stage
%
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function [BM MF Lp] = CaspPreProcInit(BM, MF, fs, fsIR);

%% calculates center frequencies of basilar-membrane filterbank in ERB
% [BM.NrChannels, BM.CenterFreq] = getGFBCenterERBs(BM.MinCF, BM.MaxCF, BM.Align, BM.Dens);
% BM.CenterFreqs = erbtofreq(BM.CenterFreq); % convert from ERBs to freqs
warning('previous line by-passed by AO and replaced by following line...')

% find the center frequencies used in the filterbank, 1 ERB spacing
BM.CenterFreqs  = erbspacebw(BM.MinCF, BM.MaxCF, BM.Dens, []);
BM.CenterFreq   = freq2erb(BM.CenterFreqs);
BM.NrChannels = length(BM.CenterFreqs);

%% calculate filter coefficients for basilar-membrane filterbank and lowpass
[BM.b(1,:,:), BM.a(1,:,:)]=getGFBFilterCoefs(BM.NrChannels, BM.CenterFreqs, BM.BW, fs);
[Lp.b1, Lp.a1] = butter(1, 1000*2/fs);

%% lowest and higest CFs of the MFB as function of Gtf.CF
MF.low = (BM.CenterFreq .* 0 + 1) .* MF.MinCF;				% set lowest mf as constant value
MF.up = min(BM.CenterFreqs .* MF.gtfac, MF.MaxCF); 	% set highest mf as proportion of gt_CF's

%% calculate coefficients and center frequencies of modulation filters
switch MF.style
    case 'mfbtd_drnl'
        [MF.CenterFreq, out] = mfbtd_drnl(1,min(MF.low),max(MF.up),MF.lpco,fsIR);
        
    case 'lp'
        [Lp.b(1,:,:), Lp.a(1,:,:)]=IRIfolp(1000,fs);
        MF.CenterFreq = 1;                     % just to provide an adequate vlaue for size-operation
        [MF.MLpb MF.MLpa] = IRIfolp(1/(2*pi*0.02)/fsIR);        % provide the modulation lp coeffs for 20ms time constant
    otherwise
        error('CaspPreProcInit: illegal MF.style')
end
MF.NrChannels = size(MF.CenterFreq,2); % calculate maximum number of modulation filters
% eof


