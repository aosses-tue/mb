function Vocoder_CISIM
% function Vocoder_CISIM
%
% 1. Description:
%       Cochlear Implant simulation (CISIM) audio files were generated 
%       manually from KUL_Sim_F0m.mdl Simulink model, by loading UW files
%       into the audio block available at the model. Then, the reconstructed
%       signal (8 channel vocoder, Gaussian noise) was stored as a MAT-variable
%       including also p, q structs (containing CI info, and set-up) and this
%       script takes the available MAT files to convert them into WAV files
%       sampled at 44.1 kHz to be used in APEX. Romain Peeters' clinical map
%       was loaded into the simulink model.
% 
%       This audio files/APEX experiment were used for Brownbag meeting 
%       celebrated on Monday 16th at HTI, TU/e
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 12/6/2014
% Last update: 12/6/2014 % Update this date manually
% Last used: 12/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputpath = 'C:\Users\aosses\Dropbox\TUe\vocoder\';
outputdir = Get_TUe_paths('outputs');

files = dir([inputpath 'UW*.mat' ]);
SPL = 60;

for i = 1:length(files)
    wavfilename = [outputdir Delete_extension(files(i).name,'mat')];
    load([inputpath files(i).name])
    voc = setdbspl(voc_out,SPL);
    voc44100 = resample(voc,44100,round(p.CFG.Fs));
    Wavwrite(voc,p.CFG.Fs,wavfilename)
    Wavwrite(voc44100,44100,[wavfilename '_Hz_CISIM']);
    try
        Wavwrite(dir_out,p.CFG.Fs,[wavfilename '-dir']);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])