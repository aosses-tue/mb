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

% TODO:
%
% Increase FrameLength to 2 seconds
% Correct displayed values at the end (see 'output').

close all

pathaudio = Get_TUe_paths('outputs');

filename = [pathaudio 'ref_fluct'];

% filename = [pathaudio 'ref_rough'];

try
    
    [x Fs] = Wavread(filename);
    
catch % if reference does not exist, then it is created...
    
    N_window = 8192; % 8192;
    dur = N_window*10;
    
    opts.bDoFluct = 1;
    opts.bDoRough = 0;
    opts.dur = (dur/44100); %  approx. 200e-3;
    opts.bGen_test_tones = 1;
    opts.bForPsySound = 1;
    
    opts.bDoRamp = 0;
    opts.dur_ramp_ms = 10;
    opts.bDoZeroPadding = 1;
    opts.dur_zero_samples = round(N_window/2); % 4096 + 4096
    
    outs = Generate_reference_sounds(opts);
end

% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_002Hz'];
filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_004Hz'];
% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_008Hz'];
% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_016Hz'];

[x Fs] = Wavread(filename);
%
% 200e-3            vacil       SPLout
% -------------------------------------
% ref_fluct         1.4295      82.4621
% ..fmod_002Hz      0.2547      73.6141
% ..fmod_004Hz      0.8201      68.2238
% ..fmod_008Hz      0.8704      73.4461
% ..fmod_016Hz      0.0271      73.3777

% 2                 vacil       SPLout
% -------------------------------------
% ref_fluct          
% ..fmod_002Hz       
% ..fmod_004Hz       
% ..fmod_008Hz       
% ..fmod_016Hz       

t = (1:length(x))/Fs;
figure;
plot(t,x)
opts = []; % we clear opts

% opts.nAnalyser = 15; % Roughness
opts.nAnalyser = 20; % Fluctuation strength
opts.CalMethod = 1; % 0 dBFS = 100 dB
% opts.Author = 'DW';

filename = [filename '.wav'];
PsySoundCL(filename,opts);

x =  [];
y =  x;

% Output:
% Analyser Number: 20
% Param 1 out of 3: Roughness
% 	 [average] = [NaN]
% 	 [min max] = [0.000 1.249]
% Analyser Number: 20
% Param 2 out of 3: Specific Roughness
% 	 [average] = [NaN]
% 	 [min max] = [NaN NaN]
% Analyser Number: 20
% Param 3 out of 3: Average Roughness
% 	 [average] = [NaN]
% 	 [min max] = [-0.024 0.149]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
