function Generate_Fastl2007_textGrid(track_no)
% function Generate_Fastl2007_textGrid(track_no)
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
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures2generate = {    9, 3,'track_34'; ...
                       10, 1,'track_35'; ...
                       11, 1,'track_38'; ...
                       11, 3,'track_39'};
                       
if nargin == 0
    
    for i = 1:size(figures2generate,1)
        tmpinfo = Check_Fastl2007(figures2generate{i,1}, figures2generate{i,2});
        Generate_Praat_textGrid([Get_TUe_paths('db_fastl2007_src') figures2generate{i,3} '.wav'],tmpinfo.timesep);

        % move to Fastl2007's directory
    end
    
else
    
    cont = length(figures2generate);
    while cont >= 1
        if strcmp(num2str(track_no),figures2generate{cont,3}(end-1:end))
            i = cont;
        end
        cont = cont - 1;
    end
    
    try
        tmpinfo = Check_Fastl2007(figures2generate{i,1}, figures2generate{i,2});
        Generate_Praat_textGrid([Get_TUe_paths('db_fastl2007_src') figures2generate{i,3} '.wav'],tmpinfo.timesep);
    catch
        error('Requested Fastl''s figure not processed yet, try another track_no or run this function without arguments');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
