function y = r20141107_fluctuation
% function y = r20141107_fluctuation
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/11/2014
% Last update on: 05/11/2014 % Update this date manually
% Last use on   : 05/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

pathaudio = Get_TUe_paths('outputs');

% filename = [pathaudio 'ref_fluct'];
filename = [pathaudio 'ref_rough'];

try
    
    Wavread(filename);
    
catch % if reference does not exist, then it is created...
    
    opts.bDoFluct = 1;
    opts.dur = (8192/44100); %  approx. 200e-3;
    opts.bGen_test_tones = 0;
    opts.bForPsySound = 1;
    
    opts.bDoRamp = 1;
    opts.dur_ramp_ms = 25;
    opts.bDoZeroPadding = 1;
    opts.dur_zero_samples = 4096; % 4096 + 4096
    
    outs = Generate_reference_sounds(opts);
end

opts = []; % we clear opts

% opts.nAnalyser = 15; % Roughness
opts.nAnalyser = 20; % Roughness
opts.CalMethod = 1; % 0 dBFS = 100 dB
% opts.Author = 'DW';

filename = [filename '.wav'];
PsySoundCL(filename,opts);

x =  [];
y =  x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
