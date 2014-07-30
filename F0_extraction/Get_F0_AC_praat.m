function [t F0] = Get_F0_AC_praat(inputfile,outputfile, type)
% function [t F0] = Get_F0_AC_praat(inputfile,outputfile, type)
%
% 1. Description:
% 
% % Example (Linux):
%   inputfile = '~/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu033_nodate/Audio_files-in-noise/Sweep-log-44100-up-to-1kHz-20dBFS-dir_out-SNR-99.wav';
%   outputfile = '/home/alejandro/Bureaublad/F0-test.txt';
%   Get_F0_AC_praat(inputfile,outputfile);
%
% % 2. Example + plot:
% 
% inputfile  = '~/Documenten/Databases/dir01-Instruments/Voice-of-dragon/03-Wav-files-predicted/03-Wav-files-calibrated/modus-4-v_2filt.wav';
% outputfile = '~/Bureaublad/F0-test.txt';
% type = 3;
% Get_F0_AC_praat(inputfile,outputfile,type);
% info.F0max = 1400; 
% [t F0] = Get_F0_praat_from_txt(outputfile, info);
% plot(t,F0)
% 
%
% Programmed by Alejandro Osses, ExpORL 2014
% Last updated on: 19/05/2014
% Last used on   : 19/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    type = 0; % As applied by Alejandro Osses, ExpORL, Physical validation
end

switch type
    
    case 0 % Alejandro Osses, ExpORL
        
        framelen    = 10e-3; % speechFrameLength
        minf0       = 75; % pitchFloor
        numcand     = 15; % numCandidates
        % veryAccurate = no % Set as constant
        silence_thr = 0.01; % silenceThreshold
        voicing_thr = 0.45; % voicingThreshold
        octave_cost = 0.01; % octaveCost
        octave_jump_cost = 0.35; % octaveJumpCost
        vUv_cost    = 0.14; % vUvCost 
        maxf0       = 400; % pitchCeil (changed from 1000 to 400 Hz on Februray 21st, 2014)
        % veryAccurate = no % Set as constant
        
    case 1 % Method 'ac' by de Cheveigné 2012
        
        framelen    = 10e-3; % speechFrameLength
        minf0       = 40; % pitchFloor
        numcand     = 15; % numCandidates
        % veryAccurate = no % Set as constant
        silence_thr = 0.0; % silenceThreshold, modified
        voicing_thr = 0.0; % voicingThreshold
        octave_cost = 0.01; % octaveCost
        octave_jump_cost = 0.0; % octaveJumpCost
        vUv_cost    = 0.0; % vUvCost 
        maxf0       = 800; 
        
    case 2
        
        framelen    = 10e-3; % speechFrameLength
        minf0       = 40; % pitchFloor
        numcand     = 15; % numCandidates
        % veryAccurate = no % Set as constant
        silence_thr = 0.01; % silenceThreshold, modified
        voicing_thr = 0.0; % voicingThreshold
        octave_cost = 0.0; % octaveCost
        octave_jump_cost = 0.0; % octaveJumpCost
        vUv_cost    = 0.05; % vUvCost 
        maxf0       = 800; 

    case 3 % Alejandro Osses, TU/e
        
        framelen    = 10e-3; % speechFrameLength
        minf0       = 75; % pitchFloor
        numcand     = 15; % numCandidates
        % veryAccurate = no % Set as constant
        silence_thr = 0.01; % silenceThreshold
        voicing_thr = 0.45; % voicingThreshold
        octave_cost = 0.01; % octaveCost
        octave_jump_cost = 0.35; % octaveJumpCost
        vUv_cost    = 0.14; % vUvCost 
        maxf0       = 1400;
    
end

if nargin < 3
    
    txt2disp = ['To pitch (ac)... '];
    txt2disp = [txt2disp num2str(framelen) ' ' num2str(minf0) ' ' num2str(numcand) ' no '...
                num2str([silence_thr voicing_thr octave_cost octave_jump_cost vUv_cost maxf0])]
            
    disp(txt2disp)
    
end

local_praat     = Get_TUe_paths('praat');
local_praat_sc  = Get_TUe_paths('praat_scripts');

if isunix
    script  = [local_praat_sc 'Get_F0_ACF_unix.praat'];
else
    script  = [local_praat_sc 'GetF0ACFwin.praat'];
end

if isunix
    command4system = [local_praat ' ' script ' "' inputfile '" "' outputfile '" "' num2str(framelen) '" "' num2str(minf0) '" "' num2str(numcand) '" "no" "' num2str(silence_thr) '" "' num2str(voicing_thr) '" "' num2str(octave_cost) '" "' num2str(octave_jump_cost) '" "' num2str(vUv_cost) '" "' num2str(maxf0) '"'];
else
    % command4system = [local_praat ' ' script]; % this command is used only to open Praat (run it then manually: Ctrl+R)
    command4system = [local_praat ' ' script ' "' inputfile '" "' outputfile '" "' num2str(framelen) '" "' num2str(minf0) '" "' num2str(numcand) '" "no" "' num2str(silence_thr) '" "' num2str(voicing_thr) '" "' num2str(octave_cost) '" "' num2str(octave_jump_cost) '" "' num2str(vUv_cost) '" "' num2str(maxf0) '"'];
end

disp([mfilename '.m: ' command4system])
[s r] = system( command4system );

if nargout ~= 0
    [t F0] = Get_F0_praat_from_txt(outputfile); % F0max constrained to 1000 Hz
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
