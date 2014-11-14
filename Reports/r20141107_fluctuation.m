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
% Last use on   : 14/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
%
% Increase FrameLength to 2 seconds
% Correct displayed values at the end (see 'output').

close all

pathaudio = Get_TUe_paths('outputs');

filename = [pathaudio 'ref_fluct'];

% filename = [pathaudio 'ref_rough'];

N = 8192;
bDebug = 1;

try
    
    [x Fs] = Wavread(filename);
    
catch % if reference does not exist, then it is created...
    
    N_window = N; % 8192;
    dur = N_window*2;
    
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

% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_002Hz_60_dBSPL'];
filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_004Hz_60_dBSPL'];
% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_008Hz_60_dBSPL'];
% filename = [pathaudio 'test_fluct_fc_1000_AM_m_100_fmod_016Hz_60_dBSPL'];
% filename = [pathaudio 'test_rough_fc_0500_AM_m_100_fmod_130Hz_60_dBSPL.wav']; % it should be 0
% filename = [pathaudio 'test_fluct_fc_1000_FM_dev_700_fmod_002Hz_70_dBSPL']
% filename = [pathaudio 'test_fluct_fc_1000_FM_dev_700_fmod_004Hz_70_dBSPL']
% filename = [pathaudio 'test_fluct_fc_1000_FM_dev_700_fmod_008Hz_70_dBSPL']
% filename = [pathaudio 'test_fluct_fc_1000_FM_dev_700_fmod_016Hz_70_dBSPL']

[x Fs] = Wavread(filename);

t = (1:length(x))/Fs;
% figure;
% plot(t,x)
opts = []; % we clear opts

% opts.nAnalyser = 15; % Roughness
opts.nAnalyser = 20; % Fluctuation strength
opts.CalMethod = 1; % 0 dBFS = 100 dB

starti = N+1;
% insig = x(starti:starti + N-1);
insig = x(starti:starti + N-1);
t = ( 1:length(insig) )/Fs;

if bDebug
    figure(1)
    plot(t,insig);
    xlabel('Time [s]')
    ylabel('Amplitude')
    title(sprintf('Test signal, level = %.2f [dB]',rmsdb(insig)+90))
end

out = FluctuationStrength_offline(insig,Fs,N,bDebug);

% filename = [filename '.wav'];
% PsySoundCL(filename,opts);

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
