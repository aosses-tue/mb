function Generate_piano_textGrid(track_no)
% function Generate_piano_textGrid(track_no)
%
% 1. Description:
%       track_no - track number
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3.1 Stand-alone example, generate all TextGrids:
%       Generate_Fastl2007_textGrid;
%
% 3.2 Stand-alone example, generate only TextGrid for track 39:
%       Generate_Fastl2007_textGrid(39);
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 14/08/2014
% Last update on: 14/08/2014
% Last use on   : 25/11/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir2lookat = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\01-Tuned-at-44100-Hz\';

registers = {   ['Dsh1' delim]; ...
                ['F1' delim]; ...
                ['C2' delim]; ...
                ['Ash2' delim]; ...
                ['F3' delim]; ...
                ['C4' delim]; ...
                ['A4' delim]; ...
                ['Csh5' delim]};
                % ['C6' delim]
                % ['G6' delim]

   
for i = 1:length(registers)
    dirtmp = [dir2lookat registers{i}];
    filenames = Get_filenames(dirtmp,'*.wav');
    for j = 1:length(filenames)
        timesep = Check_piano_sounds(filenames{j});
        Generate_Praat_textGrid([dirtmp filenames{j}],timesep);
    end
    % move to Fastl2007's directory
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
