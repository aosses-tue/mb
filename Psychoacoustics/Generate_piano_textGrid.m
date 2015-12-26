function Generate_piano_textGrid
% function Generate_piano_textGrid
%
% 1. Description:
%       track_no - track number
% 
% 2. Stand-alone example, generate all TextGrids:
%       Generate_piano_textGrid;
%  
% 3. Additional info:
%       Tested cross-platform: No
%       See also Check_piano_sounds, Generate_Praat_textGrid, Generate_Fastl2007_textGrid (similar processing)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 14/08/2015
% Last update on: 14/08/2015
% Last use on   : 22/12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir2lookat = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Tuned-at-44100-Hz' delim];

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
