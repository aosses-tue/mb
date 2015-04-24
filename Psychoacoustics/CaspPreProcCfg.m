%--------------------------------------------------------------------------
% CaspPreProcCfg.m
%--------------------------------------------------------------------------
% relevant parameters in the pre-processing for the basilar-membrane and
% modulation filterbank are assigned
%
% usage
%   [BM MF] = CaspPreProcCfg
%
% output
%   BM          : structure containing filter coefficients, center frequencies
%                 and stuff for the basilar-membrane filterbank
%   MF          : structure containing filter coefficients, center frequencies
%                 and stuff for the modulation filterbank
%
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function [BM MF] = CaspPreProcCfg

%% Basilar filterbank variables
BM.MinCF = 80;               % lowest CF in Hz
BM.MaxCF = 8000; % 10000      % highest CF in Hz
warning('Parameter changed by AO: temporarily')

BM.Align = 1000;              % base frequency in Hz
BM.BW    = 1.0;               % bandwidth of the filter in ERB
BM.Dens  = 1.0;               % filter density in 1/ERB

%% Modulation filterbank variables
MF.MinCF = 0;                 % lowest modulation filter
MF.MaxCF = 1000;              % upper most modulation filter (mj init 1000)
MF.gtfac = 0.25;         	  % factor for deriving the upper most mf depending on the gt_CF
MF.lpco = 1;                  % lowpass cut off of the modulation filter 1 = 2.5 Hz, 2 = 7.5 Hz
MF.style = 'mfbtd_drnl';      % use 'mfbtd_drnl' for multi or single channel and 'lp' for modulation lowpass only
%eof
