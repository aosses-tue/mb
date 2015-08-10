function filename = AMTControl_Examples(nExample)
% function filename = AMTControl_Examples(nExample)
%
% 1. Description:
%       Wav files are generated if no output is specified
%   
% 2. Stand-alone example:
%       % Tone generation + save to file of audios of example 1
%       AMTControl_Examples(1); 
%       
%       AMTControl_Examples(3); 
% 
%       % To get the filenames but not storing the audio files
%       filename = AMTControl_Examples(1); 

% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20150522_update (for example 1)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 10/08/2015
% Last update on: 10/08/2015 
% Last use on   : 10/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dirout = [Get_TUe_paths('outputs') 'AMTControl-examples' delim];

if nargin == 0
    nExample = 1;
end

switch nExample
    case 0 % default
        
        file1 = ['dau1996b_expI_noisemasker.wav'];
        file2 = ['dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav'];
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
    case 1 % Weber experiment
        
        fs = 44100;
        f = 1000;
        dur = 4; % in seconds
        SPLs = [60 42]; % see r20150522_update. Estimated dprime should be aroung 1.25
        
        [outsig1 file1] = Il_create_tone(f,dur,fs,SPLs(1));
        [outsig2 file2] = Il_create_tone(f,dur,fs,SPLs(2));
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
        if nargout == 0
            Wavwrite( outsig1,fs,filename{1} );
            Wavwrite( outsig2,fs,filename{2} );
        end
        
    case 2 % Fluctuation strength samples
        
        file1 = ['fluct_test_bbn_AM_m_000_fmod_004Hz_60_dBSPL.wav'];
        file2 = ['fluct_test_bbn_AM_m_070_fmod_004Hz_60_dBSPL.wav'];
        
        filename{1} = [dirout file1];
        filename{2} = [dirout file2];
        
    case 3 % Modulation detection
        
        fs = 44100;
        BW = 3;
        fc = 5000;
        fmod = 20;
        dur = 1;
        Mdept = 0.5; % well above threshold
        SPL = 70;
        Fs = 44100;
        [y env] = AM_random_noise_BW(fc,BW,SPL,dur,fs,fmod,Mdept);

        % If you want to store the output (Wav file):
        AM_random_noise_BW(fc,BW,SPL,dur,Fs,fmod,Mdept);
        
%         file1 = ['fluct_test_bbn_AM_m_000_fmod_004Hz_60_dBSPL.wav'];
%         file2 = ['fluct_test_bbn_AM_m_070_fmod_004Hz_60_dBSPL.wav'];
%         
%         filename{1} = [dirout file1];
%         filename{2} = [dirout file2];
        
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function [outsig filename] = Il_create_tone(f,dur,fs,lvl)

% Creating test tone:
y       = .5*Create_sin(f,dur,fs,0);
rampup  = 5; % ms
rampdn  = 5; % ms
insig    = Do_cos_ramp(y,fs,rampup,rampdn);

outsigtmp   = setdbspl(insig,lvl);

outsig      = [Gen_silence(200e-3,fs); ...
               outsigtmp;
               Gen_silence(200e-3,fs)];
           
filename = sprintf('tone-f-%.0f-Hz-at-%.0f-dB-dur-%.0f-s.wav',f,lvl,dur);
