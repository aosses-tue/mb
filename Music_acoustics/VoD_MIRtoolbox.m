function VoD_MIRtoolbox(filename,fs)
% function MIRtoolbox(filename,fs)
%
% 1. Description:
%       fs required if 'filename' is numeric.
%       MIR toolbox: interesting but discarded from our analysis according
%       to meeting with Armin on 23/07/2014 (no psychoacoustics). However
%       there are some very interesting tools.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       VoD_MIRtoolbox;
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/7/2014
% Last update on: 21/7/2014 % Update this date manually
% Last used on  : 21/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    
    close all
    m = Get_TUe_subpaths('db_voice_of_dragon');
    rootfolder = m.dir_calibrated_p;

    filename = [rootfolder  'modus-1-v_2filt.wav'];

end


%% Section 1: 

% Extracting a part, then playing the excerpt:
if ischar(filename)
    y = miraudio(filename, 'Extract',0,1); % extracting from second 1 to 2
    miraudio(filename, 'Extract',0,1)
else
    if isnumeric(filename)
        if nargin < 2
            warning('considering sampling frequency as 44100 Hz')
            fs = 44100;
        end
        miraudio(filename,fs,'Extract',0,1)
        y = miraudio(filename,fs,'Extract',0,1);
    end
end

%% 2. Spectral anamysis:
% Spectral decomposition

% temporal evolution, with Freq max of 3000 Hz:
mirspectrum(y,'Frame','Max',3000)


% Spectral flux:
s = mirspectrum(y,'Frame','Max',3000);
mirflux(s)
 
% Brightness:
mirbrightness(y,'Frame')

% Spectral centroid:
mircentroid(y,'Frame')
 
% % Roughness:
mirroughness(y)
% mircentroid(a,'Frame')
 
%% 4. Transcription:
% MIR onsets: detection of succesive notes
mironsets(y)
 
% Attack:
mironsets(y,'Attacks')

% % Release:
% a = miraudio('ragtime','Extract',1,2);
% % mironsets(a,'Releases') % Failing
% 
% % Attack of notes:
% mirattackslope(a)
% 
% % Density of notes
% mireventdensity('ragtime','Frame') % mireventdensity('ragtime')
% 
% % Pitch height
% mirpitch('ragtime','Frame')
% 
%% 6. Structure:
% Similarity matrix
mirsimatrix(y)
% mirsimatrix('ragtime','Frame',1) % frames of 1 second
% 
% % MFCC:
% mirmfcc('ragtime','Frame')
% 
% % MFCC, temporal evolutions in timbral domain:
% c = mirmfcc('ragtime','Frame');
% mirsimatrix(c)
% 
% Novelty curve
mirnovelty(y)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end