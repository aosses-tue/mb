
function quick_check_trial_levels
% function quick_check_trial_levels
%
% Run this script to generate the loudest and softest trials for Pitch 
% Ranking (PR), Melodic Contour Identification (MCI). LIST and Lilliput are 
% also taken into account but just when running the APEX experiment:
%   ../Meas/Experiments/Tests_XML/Testen-en-Training/TEST_Protocol-check-MAP.xml
%
% In that experiment, the trials PR_p00.wav, PR_m10.wav, MC_R4_p00.wav, 
% MC_R4_m04.wav are loaded in addition to 'wdz6.wav' (LIST) and 
% 'lach_Febr2012.wav'. Make sure that you added manually these two last trials 
% from Gilbert.
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List_of_files = {   'generate_silence.m'};
                
[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputdirname = 'Protocol-check-MAP';

try
    
    signal = [];

    % 0 dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_104_Hz.wav');
    dur_silence = 30e-3;
    s = generate_silence(Fs,dur_silence);

    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_131_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_147_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_165_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_208_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_262_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];

    [x, Fs] = wavread('PR_Stimuli/UW_LB_ACE_294_Hz.wav');
    signal = [signal; x];
    signal = [signal; s];
    
    try
        wavwrite(signal,Fs,[outputdirname '/PR_p00'])
    catch
        mkdir(outputdirname)
        wavwrite(signal,Fs,[outputdirname '/PR_p00'])
    end
    
    % -10 dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wavwrite(From_dB(-10)*signal,Fs,[outputdirname '/PR_m10'])
    disp(['PR wav file lasts : ' num2str(length(signal)/Fs)])

    % MC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    copyfile('MC_Stimuli/UW_LB_ACE_R4_104_Hz.wav',[outputdirname '/MC_R4_p00.wav']);

    signal = wavread('MC_Stimuli/UW_LB_ACE_R4_104_Hz.wav');

    wavwrite(From_dB(-4)*signal,Fs,[outputdirname '/MC_R4_m04'])
    disp(['MC wav file lasts : ' num2str(length(signal)/Fs)])
    
catch
    error([mfilename '.m: Point your MATLAB path to ../Meas/Music/'])
end

end