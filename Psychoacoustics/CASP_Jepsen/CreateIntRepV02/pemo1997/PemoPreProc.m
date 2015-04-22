%--------------------------------------------------------------------------
% PemoPreProc.m
%--------------------------------------------------------------------------
%   runs pre-processing of Pemo and calculates the internal representation
%   of the input signal with pemo
%
% usage
%   PemoPreProc(x, fs)
%
% input
%   x           : input signal in wav-format
%   fs          : sample rate of input signal
%
% output
%   IntRep      : structure containing internal representation
%   BM          : basilar filterbank coefficients & stuff
%   Mf          : modulation filterbank coefficients & stuff
%
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function [IntRep BM MF] = PemoPreProc(x, fs)

%% if no input arguments are set, the demo input signal is loaded
if nargin < 2, fs = 16000; end
if nargin < 1, x = wavread('x.wav'); end

% ------------------------------------------------------------------------
%% pemo pre-processing
% ------------------------------------------------------------------------

%% implementation stuff
% signal format one channel (i.e. one column) only
[row col] = size(x);
if col > 1
    warning('signal format: only one channel (first column) will be processed')
    x = x(:,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT TOUCH!!                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doing up- or downsampling because outer-middle ear filter coeff's are
% calculated for input signals of 44.1 kHz only in casp, done here too, to
% maintain direct comparability
if fs < 44100;
    x = resample(x,44100, fs);
    fs = 44100;
    disp('signal is upsampled to 44.1 kHz')
elseif fs > 44100
    xtmp = resample(x,fs,44100);
    IntRep.resampleLen = floor(length(x) / (44100/fs));
    fs = 44100;
    x = xtmp(1:IntRep.resampleLen);
    disp('signal is downsampled to 44.1 kHz')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% avoid truncation errors due to resampling
IntRep.resampleFac = 4;
IntRep.resampleLen = floor(length(x) / IntRep.resampleFac);
IntRep.fs = fs / IntRep.resampleFac;

% necessary parameters for basilar-membrane and modulation filterbank
[BM MF] = PemoPreProcCfg;
[BM MF Lp] = PemoPreProcInit(BM, MF, fs, IntRep.fs);

% provide an adequate structure
IntRep.x = x;
IntRep.Y = zeros(IntRep.resampleLen,BM.NrChannels,MF.NrChannels);

%% calculate internal representation
for ChannelNr = 1:BM.NrChannels		% gammatone filter loop
    
    % stage1: gammatone filterbank
    bc = squeeze(BM.b(1,ChannelNr,:));
    ac = squeeze(BM.a(1,ChannelNr,:));
    y = 2*real(filter(bc,ac,x));
    IntRep.BM(:,ChannelNr) = y;
    
    % stage2: haircell model
    y = max(y,0); %hwr
    y = filter(Lp.b1,Lp.a1,y);
    IntRep.hc(:,ChannelNr) = y;
    
    % stage3: adaptation loops
    y = max(y,1e-5);
    y = nlal_lim(y,fs,0,1e-5);
    IntRep.adapt(:,ChannelNr) = y;
    
    % downsampling
    y = resample(y,1,IntRep.resampleFac);
    y = y(1:IntRep.resampleLen);
    
    % stage4: modulation filterbank
    switch MF.style
        case 'mfbtd'
            [MF.CenterFreq, y] = mfbtd(y,MF.low(ChannelNr),MF.up(ChannelNr),MF.lpco,IntRep.fs); % calculate coefficients and center frequencies of modulation filters
            y = mfbtdpp(y, MF.CenterFreq,IntRep.fs);	% postprocessing for mfb
        case 'lp'
            MF.CenterFreq = 1;
            y = filter(MF.MLpb, MF.MLpa, y);
        otherwise
            error('PemoPreProc: illegal MF.style')
    end
    
    % adds actual channel to the internal representation matrix
    IntRep.Y(:,ChannelNr,1:length(MF.CenterFreq)) = y;
end

%eof


