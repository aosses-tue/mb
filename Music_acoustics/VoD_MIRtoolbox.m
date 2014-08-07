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
% 3.1 Stand-alone example:
%       % Make sure you have not run PsySound3
%       % Make sure MIRToolbox v. 1.5 has been added to path
%       VoD_MIRtoolbox;
%
% % Audio excerpt into double:
%       xx = get(y,'Data');
%       xx = xx{1};
%       xx = xx{1};
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/07/2014
% Last update on: 06/08/2014 % Update this date manually
% Last use on   : 06/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    
    close all
    m = Get_TUe_subpaths('db_voice_of_dragon');
    rootfolder = m.dir_calibrated_m;

    filename = [rootfolder  'modus-1_v2-2filt.wav']; % Ac.mode = 2, measured signal
    
    disp('In case is not running, please remove MIR tb from PsySound and add MIR toolbox v15')
    disp('press any button to continue...')
    pause
    % addpath('D:\MATLAB_git\tb_MIR_v15\MIRToolbox\')
end


misc    = Get_VoD_params(0);
ac_mode = 2;
mode_idx= ac_mode - 1;
ti      = misc.ti_measured(mode_idx);
nPeriods = 3;
tf      = ti + nPeriods*misc.Tmodel(mode_idx);
%% Section 1: 
% Extracting a part, then playing the excerpt:
if ischar(filename)
    % y = miraudio(filename, 'Extract',0,1); % extracting from second 1 to 2
    y = miraudio(filename, 'Extract',ti,tf); % extracting from second 1 to 2
    mirplay(y);
    
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

mirfluctuation(y)
mirfluctuation(y,'Summary')

%% 6. Structure:
% Similarity matrix
mirsimatrix(y)
% mirsimatrix('ragtime','Frame',1) % frames of 1 second

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