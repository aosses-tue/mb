function timesep = Check_piano_sounds(filename)
% function timesep = Check_piano_sounds(filename)
%
% 1. Description:
%       This script os used to generate piano textGrid (to be used in the 
%       Praat analyser.
% 
% 2. Stand-alone example:
%       filename = 'GRAF28-Dd1.wav';
%       Check_piano_sounds(filename);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 24/11/2015
% Last update on: 25/11/2015 
% Last use on   : 29/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename = 'track_04';
% timesep = [ 1 5 7 11 13 17 39 43 45 49 51 55]; % sound between 1-5, 7-11

% paths = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\01-Tuned-at-44100-Hz\';
%% F1:
switch filename
    
    %%% Dsh1 
    case 'GRAF28-Dd1.wav'
        timesep = [0.05 2 2.05 4.5 4.55 7.3];
    case 'JBS36-Dd1.wav'
        timesep = [0 7];
    case 'JBS50-Dd1.wav'
        timesep = [0.01 2.6 3 6.5 7.3 11];
    case 'JBS50-Dd1-V2-p.wav'
        timesep = [0.01 3.2 4.4 6.3];
    case 'JBS51-4486-Dd1.wav'
        timesep = [0.03 6.8];
    case 'JBS51-4544-Dd1.wav'
        timesep = [0.04 5.5 7.3 13.8 14.7 19.9 21.5 26.7];
    case 'JBS73-Dd1.wav'
        timesep = [0.1 5];

    %%% F1:
    case 'GH05-F1.wav'
        timesep = [ 0.3 9];
    case 'GRAF28-F1.wav'
        timesep = [ 0 1.9 1.95 4.1 4.15 6.5]; 
    case 'JBS36-F1.wav'
        timesep = [ 0.1 7.8]; % longer than expected 
    case 'JBS51-4486-F1.wav'
        timesep = [0.05 8.5]; % longer than expected
    case 'NS19-F1.wav'
        timesep = [0.05 4.85 4.9 9.6];
    
    %%% C2:
    case 'GH05-C2.wav'
        timesep = [0.15 6];
    case 'GRAF28-C2.wav'
        timesep = [0.27 2.75 2.8 5.28 5.3 8.5];
    case 'JBS36-C2.wav'
        timesep = [0.01 3.45 3.5 8.25];
    case 'JBS50-C2.wav'
        timesep = [0.01 4 4.52 9.4 10.2 15];
    case 'JBS51-4486-C2.wav'
        timesep = [0.01 6.7];
    case 'JBS51-4544-C2.wav'
        timesep = [0.01 4 5 9.6 10.75 16 16.65 20.15];
    case 'JBS73-C2.wav'
        timesep = [0.01 5.01];
    case 'JBS73-C2-V2.wav'
        timesep = [0.25 5.3 5.4 8.8 9 12 12.2 15.35];
    case 'NS19-C2.wav'
        timesep = [0.01 2.6 2.7 5.15];
        
    %%% Ash2:
    case 'GH05-Ad2.wav'
        timesep = [0.41 6];
    case 'GRAF28-Ad2.wav'
        timesep = [0.01 1.85 1.9 4 4.1 6.2];
    case 'JBS36-Ad2.wav'
        timesep = [0.015 5 5.4 9.8];
    case 'JBS50-Ad2.wav'
        timesep = [0.01 3.3 3.4 7.3 7.6 11 12 15.5];
    case 'JBS50-Ad2-V2-p.wav'
        timesep = [0.15 1.7];
    case 'JBS51-4486-Ad2.wav'
        timesep = [0.03 6.5];
    case 'JBS51-4544-Ad2.wav'
        timesep = [0.01 3.3 3.4 7.3 8 11.5];
    case 'JBS73-Ad2.wav'
        timesep = [0.05 4.5];
    case 'NS19-Ad2.wav'
        timesep = [0.1 2.75 2.8 5.03];
            
    %%% F3:    
    case 'GH05-F3.wav'
        timesep = [0.07 4.15 4.2 9.15];
    case 'GRAF28-F3.wav'
        timesep = [0.01 1.85 3.85 3.9 5.85];
    case 'JBS36-F3.wav'
        timesep = [0.1 3.4 3.5 6.5];
    case 'JBS50-F3.wav'
        timesep = [0.01 3.8 3.9 7.3 7.6 11.5 11.55 14.85];
    case 'JBS51-4486-F3.wav'
        timesep = [0.01 3.25 3.3 6.5];
    case 'JBS51-4544-F3.wav'
        timesep = [0.01 3.1 3.2 7 7.05 11.3 11.7 16.5];
    case 'JBS73-F3.wav'
        timesep = [0.01 2.7];
    case 'NS19-F3.wav'
        timesep = [0.01 1.85 1.9 3.8];
               
    %%% C4:    
    case 'GH05-C4.wav'
        timesep = [0.15 3.3 3.4 8];
    case 'GRAF28-C4.wav'
        timesep = [0.15 3.2 3.3 5.5 6.2 8.5];
    case 'JBS36-C4.wav'
        timesep = [0.18 3.6 3.65 7.5];
    case 'JBS50-C4.wav'
        timesep = [0.01 3.6 4 8 8.65 12.7];
    case 'JBS51-4486-C4.wav'
        timesep = [0.05 3.25 3.3 8];
    case 'JBS51-4544-C4.wav'
        timesep = [0.01 3 3.35 7.2 7.3 10.5 11 14.2];
    case 'JBS73-C4.wav'
        timesep = [0.01 2.5 2.6 6.1 6.2 9.6 9.7 12];
    case 'NS19-C4.wav'
        timesep = [0.1 3.1 3.15 6.3 6.35 9.41 9.5 12.25 12.3 15 15.1 15.15 17.5];
        
	%%% A4:    
    case 'GH05-A4.wav'
        timesep = [0.3 3.5 3.9 7.3];
    case 'GRAF28-A4.wav'
        timesep = [0.05 1.7 1.75 3.8 3.85 5.9];
    case 'JBS36-A4.wav'
        timesep = [0.18 3.7 3.8 7.9];
    case 'JBS50-A4.wav'
        timesep = [0.01 3.3 3.35 7.5 7.8 12 12.1 16.4];
    case 'JBS51-4486-A4.wav'
        timesep = [0.03 3.7 3.8 6.8];
    case 'JBS51-4544-A4.wav'
        timesep = [0.01 3.4 3.45 7.5 7.9 11.5 12.2 15];
    case 'JBS73-A4.wav'
        timesep = [0.01 4.2];
    case 'NS19-A4.wav'
        timesep = [0.01 1.6 1.65 3.25];

	%%% Csh5:    
    case 'GH05-Cd5.wav'
        timesep = [0.5 3.5 3.6 7.3 7.5 11.5];
    case 'GRAF28-Cd5.wav'
        timesep = [0.01 2.1 2.15 4.3 4.35 6.6];
    case 'JBS36-Cd5.wav'
        timesep = [0.15 3.6 3.7 3.8 7.1 9.4];
    case 'JBS51-4486-Cd5.wav'
        timesep = [0.1 2.9 3 5.9 6 8.2];
    case 'JBS51-4544-Cd5.wav'
        timesep = [0.01 2.7 3.1 6.5 6.9 9.5];
    case 'JBS51-4544-Cd5-V2.wav'
        timesep = [0.01 3.4 3.5 7.5 7.6 11.5 11.9 14.7];
    case 'JBS73-Cd5.wav'
        timesep = [0.01 2.65];
    case 'NS19-Cd5.wav'
        timesep = [0.05 1.7 1.75 3.25];   
	case 'NS19-Cd5-V2.wav'
        timesep = [0.01 1.7 1.75 3.2]; 
        
	%%% C6
    case 'GH05-C6.wav'
        timesep = [0.2 3.5 4.1 6.6];
    case 'GRAF28-C6.wav'
        timesep = [0.01  2.1 2.15 3.6 4.15 5.9];
    case 'JBS36-C6.wav'
        timesep = [0.05 1.28 1.30 2.62 2.65 3.87];
    case 'JBS51-4486-C6.wav'
        timesep = [0.1 2.3 2.55 4.6 5.15 7.65];
    case 'JBS73-C6.wav'
        timesep = [0.06 1.92];
    case 'NS19-C6.wav'
        timesep = [0.01 1.25 1.35 2.55];   
        
    %%% G6
    
    otherwise
        disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
