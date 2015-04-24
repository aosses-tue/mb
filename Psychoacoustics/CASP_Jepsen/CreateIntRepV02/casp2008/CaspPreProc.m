%--------------------------------------------------------------------------
% CaspPreProc.m
%--------------------------------------------------------------------------
%   runs pre-processing of Pemo and calculates the internal representation
%   of the input signal with casp
%
% usage
%   CaspPreProc(x, fs)
%
% input
%   x           : input signal in wav-format
%   fs          : sample rate of input signal
%
% output
%   IntRep      : structure containing internal representation
%   BM          : basilar filterbank coefficients & stuff
%   MF          : modulation filterbank coefficients & stuff
%
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function [IntRep BM MF] = CaspPreProc(x, fs)

%% if no input arguments are set, the demo input signal is loaded
if nargin < 2, fs = 44100; end
if nargin < 1, x = wavread('x.wav'); end

% ------------------------------------------------------------------------
%% casp pre-processing
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
[BM MF]     = CaspPreProcCfg;
[BM MF Lp]  = CaspPreProcInit(BM, MF, fs, IntRep.fs);

% provide an adequate structure
IntRep.x = x;
IntRep.Y = zeros(IntRep.resampleLen,BM.NrChannels,MF.NrChannels);

%% calculates internal representation

% stage1: Outer-Middle ear filter
xStapes = OuterMiddleFilter(IntRep.x);

tmp = [];

for ChannelNr = 1:BM.NrChannels		% drnl filter loop
    % stage2: drnl filterbank
    y = drnl(xStapes',BM.CenterFreqs(ChannelNr),fs)';
    IntRep.BM(:,ChannelNr) = y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stage3: haircell model
    y = max(y,0); % HWR: half wave rectification
    y = filter(Lp.b1,Lp.a1,y);
    IntRep.hc(:,ChannelNr) = y;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % linear gain to fit ADloop operating point
    y = y*10^(50/20);
    
    % stage4: expansion
    y = y.^2;
    
    % stage5: adaptation loops
    y = max(y,1e-5);
    y = nlal_lim(y,fs,10,1e-5); % Limit 0
    warning('Overshoot changed to 10 by AO, according what it is said in Jepsen 2008')
    IntRep.adapt(:,ChannelNr) = y;
    
    tmp = [tmp y];
    
    % downsampling
    y = resample(y,1,IntRep.resampleFac);
    y = y(1:IntRep.resampleLen);
    
    % stage6: modulation filterbank
    switch MF.style
        case 'mfbtd_drnl'
            [MF.CenterFreq,y] = mfbtd_drnl(y,MF.low(ChannelNr),MF.up(ChannelNr),MF.lpco,IntRep.fs);
            y = mfbtdpp_drnl(y,MF.CenterFreq,IntRep.fs);	% postprocessing for mfb
        case 'lp'
            MF.CenterFreq = 1;
            y = filter(MF.MLpb, MF.MLpa, y);
        otherwise
            error('CaspPreProc: illegal MF.style')
    end
    % adds actual channel to the internal representation matrix
    IntRep.Y(:,ChannelNr,1:length(MF.CenterFreq)) = y;
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
