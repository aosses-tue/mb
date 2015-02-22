function [N, main_N, spec_N] = ch_dlm(sig, HL, k)
% function [N, main_N, spec_N] = ch_dlm(sig, HL, k)
%
% 1. Description:
%       Dynamic Loudness Model (Chalupper 2001): calculates loudness N, 
%       main loudness main_N and specific loudness spec_N for a signal sig 
%       and a given hearing loss (HL) (optional parameter: 
%       default 0 dB)
%       Optionally, also a k-vector can be entered (default: k=0.8) HL and k 
%       are 1x24 vectors according to Zwicker's critical bands (regarding 
%       definition of center frequencies and bandwidth)
%
%       Calibration: (fs=44.1 kHz, 107 dB SPL FS RMS(i.e. a sinusoid with 
%       amplitude of 1 has 107 dB SPL), easily 114 dB SPL = 0 dBFS
% 
% References:
% [1]   Chalupper, J.,Fastl, H. (2002): Dynamic loudness model (DLM) for 
%       normal and hearing-impaired listeners. Acustica, 88: 378-386
% [2]   Chalupper, J. (2001) - in german - : Perzeptive Folgen von
%       Innenohrschwerhoerigkeit: Modellierung, Simulation und 
%       Rehabilitation. Dissertation at the Technical University of Munich, 
%       Shaker Verlag.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000
% Edited on     : 06/01/2007 (new version with comments and examples)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 21/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Altered by MFFM Matt Flax <flatmax> for the Psy-Sound project
% Jan. 2007
% Comments :
% This file is altered from its original form to allow block based processing. 
% Block processing is a requirement of the Psy-Sound project.
% Block based processing requires that all necessary filter and data states
% are remembered between 'dlm' function calls. In order for this to be
% possible, the 'fileHandle' structure is altered accordingly.

f_abt = 1/2e-3; % 2 ms sampling period
if nargin == 2 % return the window size
  [t_pa,w,t_sb,t_sa,t_pb] = ch_staticParamDLM.m;
  [h,t, erd] = ch_tep_window(t_pb,t_pa,t_sb,t_sa,w,fs);
  N = length(h);
  
  out = N;
  return
end

if nargin < 2
    HL = zeros(1,24);
end

if nargin < 3
    k  = 0.8;
end

fs = 44100;

% Splitting of hearing loss into:
%       - HL_ihc: inner ear hair cell hearing loss
%       - HL_ohc: nonlinear component of hearing loss
HL_ohc = k.*HL;
HL_ihc = HL-HL_ohc;

% Approximation of the transfer function through the human outer and middle ear
% HPF with cut-off at (around) 65 Hz
[b, a] = ch_butter_hp(fs); % generate the butterworth filter coeffs

% filter state vector.
Z = [];

% Calculation of coefficients of critical band filterbank
S = ch_make_fttbank1(fs);

kern_l = [];

% Smoothed critical band loudness fitler creation
[smooth.b, smooth.a] = ch_int_tp(f_abt);
smooth.Zfa = [];
smooth.Zfb = [];

[N, main_N, spec_N] = run(sig);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN nested function
    function [N, main_N, spec_N] = run(sig)
  
    % Run the butterworth filter, HPF, fc = 65 Hz
    [sig, Z] = filter(b, a, sig, Z);
  
    % Applying critical band filterbank
    [fgrp, S] = ch_ftt_bank1(sig, S, f_abt,fs);
    fgrp_d    = ch_damp_a0(fgrp, HL_ihc); % Attenuation due to outer & middle
                                          % ear and inner hair cell hearing
                                          % loss
  
    % Calculation of main loudness
    kern_l = [kern_l; ch_kernlaut24_two(fgrp_d, HL_ohc)];
    
  % ii = 10; figure; plot(kern_l(:,ii)); % to visualise post-masking 
  
  % Calculation of forward masking (aka "post masking")
  try
        kern_dyn = ch_post_maskn(kern_l, f_abt);
  catch
        % caught the case where no postprocessing is possible, use a string to
        % indicate to the calling function.
        N = 'no postprocessing';
        main_N = 0;
        spec_N = 0;
        return
  end

  kern_l = []; % On successful post-processing, clear the processed data

  % Calculation of spectral masking and spectral summation of specific loudness 
  [spec_N, lauth] = ch_flankenlautheit24(kern_dyn);
  
  % Calculation of critical band loudness
  kl = ch_bark_sum(spec_N);
  
  % Smoothed critical band loudness
  [main_N, smooth.Zfa] = filter(smooth.b, smooth.a, kl, smooth.Zfa);
  main_N(find(main_N <0)) = 0;

  % Loudness integration
  [N, smooth.Zfb] = filter(smooth.b, smooth.a, lauth, smooth.Zfb);
  N(find(N <0)) = 0;
  
  end % run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 
