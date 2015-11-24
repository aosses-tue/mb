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
                ['F1' delim]};

   
for i = 1:length(registers)
    dirtmp = [dir2lookat registers{i}];
    filenames = Get_filenames(dirtmp,'*.wav');
    for j = 1:length(filenames)
        timesep = Check_piano_sounds(filenames{j});
        Generate_Praat_textGrid([dirtmp filenames{j}],timesep);
    end
    % move to Fastl2007's directory
end
    
%%%    
%     cont = length(figures2generate);
%     while cont >= 1
%         if strcmp(num2str(track_no),figures2generate{cont,3}(end-1:end))
%             i = cont;
%         end
%         cont = cont - 1;
%     end
%     
%     try
%         tmpinfo = Check_Fastl2007(figures2generate{i,1}, figures2generate{i,2});
%         Generate_Praat_textGrid([Get_TUe_paths('db_fastl2007_src') figures2generate{i,3} '.wav'],tmpinfo.timesep);
%     catch
%         error('Requested Fastl''s figure not processed yet, try another track_no or run this function without arguments');
%     end
%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
